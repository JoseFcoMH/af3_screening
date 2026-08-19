import polars as pl


def predictions_for_pair(wc) -> list[str]:
    predictions = read_pair_manifest(bait=wc.bait, prey=wc.prey)
    return expand(
        "<database>/{prediction}.tar.zst",
        prediction=predictions,
    )


rule get_pair_confidences:
    input:
        archived=predictions_for_pair,
    output:
        out_dir=temp(directory("<pair_scratch>/confidences"), group_jobs=True),
    log:
        "<logs>/get_pair_confidences/<pair_name>.log",
    priority: -50
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail
        mkdir -p {output.out_dir:q}
        printf '%s\0' {input.archived:q} \
            | xargs -0 -r -n1 -P {threads} sh -c '
                out_dir="$1"
                archive="$2"
                tar --extract --zstd \
                    --file "$archive" \
                    --directory "$out_dir" \
                    "./$(basename "$archive" .tar.zst)_confidences.json"
            ' sh {output.out_dir:q}
        """


rule get_pair_models:
    input:
        archived=predictions_for_pair,
    output:
        out_dir=temp(directory("<pair_scratch>/models"), group_jobs=True),
    log:
        "<logs>/get_pair_models/<pair_name>.log",
    priority: -50
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail
        mkdir -p {output.out_dir:q}
        printf '%s\0' {input.archived:q} \
            | xargs -0 -r -n1 -P {threads} sh -c '
                out_dir="$1"
                archive="$2"
                tar --extract --zstd \
                    --file "$archive" \
                    --directory "$out_dir" \
                    "./$(basename "$archive" .tar.zst)_model.cif"
            ' sh {output.out_dir:q}
        """


rule get_pair_summaries:
    input:
        archived=predictions_for_pair,
    output:
        out_dir=temp(directory("<pair_scratch>/summaries"), group_jobs=True),
    log:
        "<logs>/get_pair_summaries/<pair_name>.log",
    priority: -50
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail
        mkdir -p {output.out_dir:q}
        printf '%s\0' {input.archived:q} \
            | xargs -0 -r -n1 -P {threads} sh -c '
                out_dir="$1"
                archive="$2"
                tar --extract --zstd \
                    --file "$archive" \
                    --directory "$out_dir" \
                    "./$(basename "$archive" .tar.zst)_summary_confidences.json"
            ' sh {output.out_dir:q}
        """


# rule get_predictions:
#     input:
#         archived=branch(
#             lambda wc: wc.prediction in DB_PREDS,
#             [],
#             archive_batch,
#         ),
#     output:
#         confidences=temp(
#             "<pair_scratch>/{prediction}/{prediction}_confidences.json",
#             group_jobs=True,
#         ),
#         summary_confidences=temp(
#             "<pair_scratch>/{prediction}/{prediction}_summary_confidences.json",
#             group_jobs=True,
#         ),
#         model=temp(
#             "<pair_scratch>/{prediction}/{prediction}_model.cif",
#             group_jobs=True,
#         ),
#     log:
#         "<logs>/get_predictions/<pair_name>/{prediction}.log",
#     priority: -50
#     params:
#         out_dir=lambda wc, output: str(Path(output.confidences).parent),
#         database=DATABASE,
#     shell:
#         r"""
#         exec > >(tee {log:q}) 2>&1
#         set -euo pipefail
#
#         tar --extract --zstd \
#             --file {params.database:q}/{wildcards.prediction}.tar.zst \
#             --directory {params.out_dir:q} \
#             ./{wildcards.prediction}_confidences.json \
#             ./{wildcards.prediction}_model.cif \
#             ./{wildcards.prediction}_summary_confidences.json
#         """
#


rule get_and_plot_metrics:
    input:
        confs=rules.get_pair_confidences.output.out_dir,
        models=rules.get_pair_models.output.out_dir,
        summaries=rules.get_pair_summaries.output.out_dir,
    output:
        metrics_table="<pair_output>/screen_metrics.tsv",
    log:
        "<logs>/get_and_plot_metrics/<pair_name>.log",
    group:
        "stats"
    conda:
        ENVS["base"]
    params:
        ipsae_dir=lambda wc: f"{SCRATCH}/{wc.bait}__{wc.prey}/ipSAE",
        work_dir=lambda wc: f"{SCRATCH}/{wc.bait}__{wc.prey}/outputs",
        out_dir=lambda wc, output: str(Path(output.metrics_table).parent),
        pair=lambda wc: f"{wc.bait}__{wc.prey}",
        n_files=lambda wc, input: len(list(Path(input.models).glob("*.cif"))),
        models=lambda wc, input: [str(m) for m in Path(input.models).glob("*.cif")],
        script_dir=SCRIPT_DIR,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        mkdir -p {params.ipsae_dir:q} {params.work_dir:q}

        echo "Calculating metrics for {params.pair} ({params.n_files} predictions)"

        find {input.models:q} \
            -maxdepth 1 -type f -name '*.cif' -print0 \
            | xargs -0 -r -n1 -P {threads} bash -c '
                ipsae_script="$1"
                work_dir="$2"
                confs_dir="$3"
                summaries_dir="$4"
                ipsae_dir="$5"
                model="$6"

                pair_name="$(basename "$model" _model.cif)"
                pair_dir="$work_dir/$pair_name"
                mkdir -p "$pair_dir"

                model_pair="$pair_dir/${{pair_name}}_model.cif"
                jpair="$pair_dir/${{pair_name}}_confidences.json"
                summpair="$pair_dir/${{pair_name}}_summary_confidences.json"

                ln -s "$model" "$model_pair"
                ln -s "$confs_dir/${{pair_name}}_confidences.json" "$jpair"
                ln -s "$summaries_dir/${{pair_name}}_summary_confidences.json" "$summpair"

                python "$ipsae_script" \
                    "$jpair" "$model_pair" 10 10 "$ipsae_dir" \
                    || echo "WARNING: ipSAE calculation failed for $pair_name"

                rm -rf "$pair_dir"
            ' _ \
                {params.script_dir:q}/ipsae.py \
                {params.work_dir:q} \
                {input.confs:q} \
                {input.summaries:q} \
                {params.ipsae_dir:q}

        python {params.script_dir:q}/ipsae_viz_af3.py \
            --ipsae_dir {params.ipsae_dir:q} \
            --out_dir {params.out_dir:q}

        rm -rf {params.ipsae_dir:q}
        """


rule project_screen:
    input:
        confs=rules.get_pair_confidences.output.out_dir,
    output:
        out_table="<pair_output>/projected_icp.tsv",
    log:
        "<logs>/project_screen/<pair_name>.log",
    group:
        "stats"
    conda:
        ENVS["base"]
    params:
        output_dir=lambda wc, output: str(Path(output.out_table).parent),
        conf_files=lambda wc, input: [str(c) for c in Path(input.confs).glob("*.json")],
        fake_a3fdirs=lambda wc: f"{SCRATCH}/.rundata/{wc.bait}__{wc.prey}/project_screen",
        temp_manifest=lambda wc: f"{SCRATCH}/.rundata/{wc.bait}__{wc.prey}/projection_manifest.txt",
        script_dir=SCRIPT_DIR,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        fake_root={params.fake_a3fdirs:q}
        manifest={params.temp_manifest:q}

        mkdir -p "$(dirname "$manifest")" "$fake_root"
        : >"$manifest"

        trap 'rm -f "$manifest"; rm -rf "$fake_root"' EXIT

        for conf in {params.conf_files:q}; do
            name="$(basename "$conf" _confidences.json)"
            fakedir="$fake_root/$name"

            mkdir -p "$fakedir"

            conf="$(realpath "$conf")"
            link="$fakedir/$(basename "$conf")"

            ln -s "$conf" "$link"
            echo "$link" >>"$manifest"
        done

        python {params.script_dir:q}/project_screen.py \
            --predictions_list "$manifest" \
            --output_dir {params.output_dir:q} \
            --prey {wildcards.prey:q} \
            --n_cores {threads}
        """


rule summarize_all:
    input:
        tables=expand(
            f"{PAIR_OUTPUT}/screen_metrics.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
    output:
        summary_table="<work>/summary_metrics.tsv",
    log:
        "<logs>/summarize_screens.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            logger.info(f"Found {len(input.tables)} tables to concatenate.")
            tables = [pl.scan_csv(t, separator="\t") for t in input.tables]
            summary = (
                pl.concat(tables)
                .unique("pair", maintain_order=True)
                .sort("ipSAE", descending=True)
                .collect()
            )
            summary.write_csv(output.summary_table, separator="\t")
            logger.info(f"Wrote a total of {len(summary)} records.")
            logger.info("Done.")

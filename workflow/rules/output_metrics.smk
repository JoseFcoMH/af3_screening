import polars as pl


rule get_predictions:
    input:
        pred=branch(
            lambda wc: wc.prediction in DB_PREDS,
            before_update("<database>/{prediction}.tar.zst"),
            "<database>/{prediction}.tar.zst",
        ),
    output:
        confidences=temp(
            "<pair_scratch>/{prediction}/{prediction}_confidences.json",
            group_jobs=True,
        ),
        summary_confidences=temp(
            "<pair_scratch>/{prediction}/{prediction}_summary_confidences.json",
            group_jobs=True,
        ),
        model=temp(
            "<pair_scratch>/{prediction}/{prediction}_model.cif",
            group_jobs=True,
        ),
    log:
        "<logs>/get_predictions/<pair_name>/{prediction}.log",
    priority: -50
    params:
        out_dir=lambda wc, output: str(Path(output.confidences).parent),
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        tar --extract --zstd \
            --file {input.pred:q} \
            --directory {params.out_dir:q} \
            ./{wildcards.prediction}_confidences.json \
            ./{wildcards.prediction}_model.cif \
            ./{wildcards.prediction}_summary_confidences.json
        """


rule get_and_plot_metrics:
    input:
        confs=confidences_for_pair,
        models=models_for_pair,
        summaries=summaries_for_pair,
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
        out_dir=lambda wc, output: str(Path(output.metrics_table).parent),
        pair=lambda wc: f"{wc.bait}__{wc.prey}",
        n_files=lambda wc, input: len(input.models),
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        mkdir -p {params.ipsae_dir:q}

        echo "Calculating metrics for {params.pair} ({params.n_files} predictions)"

        printf '%s\0' {input.models:q} \
            | xargs -0 -r -n 1 -P {threads} bash -c '
                ipsae_script="$1"
                model="$2"

                pair_dir=$(dirname "$model")
                pair_name=$(basename "$pair_dir")
                jpair="$pair_dir/${{pair_name}}_confidences.json"

                # echo "Calculating metrics for $pair_name"
                python "$ipsae_script" "$jpair" "$model" 10 10 {params.ipsae_dir:q} || echo "WARNING: ipSAE calculation failed for $pair_name"
            ' _ "{SCRIPT_DIR}/ipsae.py"

        python {SCRIPT_DIR}/ipsae_viz_af3.py \
            --ipsae_dir {params.ipsae_dir:q} \
            --out_dir {params.out_dir:q}

        rm -rf {params.ipsae_dir:q}
        """


rule project_screen:
    input:
        confs=confidences_for_pair,
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
        prey=lambda wc: wc.prey,
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        python {SCRIPT_DIR}/project_screen.py \
            --predictions_list {input.confs} \
            --output_dir {params.output_dir} \
            --prey {params.prey} \
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
    conda:
        ENVS["base"]
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

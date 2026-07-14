from snakemake.logging import logger
import logging
from pathlib import Path
import requests
import json
import shutil
from shutil import copy, copytree

import polars as pl

from scripts.load_utils import (
    fasta_to_dict,
    fetch_uniprot_fasta,
    get_uniprot_ids_from_struc,
    dimer_template,
)

SNAKEFILE_DIR = Path(workflow.basedir).resolve()

logger.info(f"Snakefile at: {SNAKEFILE_DIR}")

SCRIPT_DIR = str(SNAKEFILE_DIR / "scripts")
ENV_DIR = str(SNAKEFILE_DIR / "envs")
CONFIG_DIR = str(SNAKEFILE_DIR / "config")


configfile: f"{CONFIG_DIR}/config.yaml"


STRUCTURE_BAITS = []
PROTEINLIST_BAITS = []

if prots := config.get("reference_proteins"):
    PROTEINLIST_BAITS = Path(prots).read_text().strip().splitlines()
if structs := config.get("reference_structures"):
    STRUCTURE_BAITS = Path(structs).read_text().strip().splitlines()

DB_PATHS = list(Path(config["predictions_database"]).glob("*.tar.zst"))
DB_PREDS = [p.stem.split(".")[0] for p in DB_PATHS]


BAITS = STRUCTURE_BAITS + PROTEINLIST_BAITS
BAITS_DICT = {Path(b).stem: b for b in BAITS}
BAITS_NAMED = list(BAITS_DICT.keys())

PREYS = Path(config.get("proteins_screened")).read_text().strip().splitlines()

WORK_DIR = config.get("work_dir")

PAIR_DIR = f"{WORK_DIR}/{{bait}}__{{prey}}"

logger.info(f"Got PREYS: {PREYS}")
logger.info(f"Got BAITS: {BAITS_NAMED}")
logger.info(
    f"Archived predictions in database (showing first 10, out of {len(DB_PATHS)}): {DB_PATHS[:10]}"
)


SCRATCH = config.get("scratch_dir", config["work_dir"] + "/scratch")
DATABASE = config["predictions_database"]


pathvars:
    database=DATABASE,
    scratch=SCRATCH,
    logs=config["work_dir"] + "/logs",


def modules_for(rule_name):
    mods = config.get("modules", {}).get(rule_name, [])
    return "" if not mods else "module load " + " ".join(mods)


def read_manifest(bait, prey):
    manifest_path = checkpoints.make_manifest.get(bait=bait, prey=prey).output.manifest
    return Path(manifest_path).read_text().strip().splitlines()


def confidences_for_pair(wc):
    predictions = read_manifest(bait=wc.bait, prey=wc.prey)
    return expand(
        "<scratch>/prediction_data/{prediction}/{prediction}_confidences.json",
        prediction=predictions,
    )


def models_for_pair(wc):
    predictions = read_manifest(bait=wc.bait, prey=wc.prey)
    return expand(
        "<scratch>/prediction_data/{prediction}/{prediction}_model.cif",
        prediction=predictions,
    )


def summaries_for_pair(wc):
    predictions = read_manifest(bait=wc.bait, prey=wc.prey)
    return expand(
        "<scratch>/prediction_data/{prediction}/{prediction}_summary_confidences.json",
        prediction=predictions,
    )


def predictions_for_pair(bait, prey):
    predictions = read_manifest(bait=bait, prey=prey)
    return expand(
        "<database>/{prediction}.tar.zst",
        prediction=predictions,
    )


def all_predictions(wc):
    preds = []
    for bait in BAITS_NAMED:
        for prey in PREYS:
            pred_names = read_manifest(bait=bait, prey=prey)
            preds.extend([f"{DATABASE}/{n}.tar.zst" for n in pred_names])
    return list(set(preds))


rule all:
    input:
        all_predictions,
        expand(
            f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
            bait=[b for b in BAITS_NAMED if BAITS_DICT[b] in STRUCTURE_BAITS],
            prey=PREYS,
        )
        if structs is not None
        else [],
        expand(
            f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.cxs",
            bait=[b for b in BAITS_NAMED if BAITS_DICT[b] in STRUCTURE_BAITS],
            prey=PREYS,
        )
        if structs is not None
        else [],
        expand(
            f"{PAIR_DIR}/processed_output/screen_metrics.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
        expand(
            f"{PAIR_DIR}/processed_output/projected_icp.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
        f"{WORK_DIR}/summary_metrics.tsv",


if PROTEINLIST_BAITS:

    rule make_baits_from_proteinlist:
        input:
            source=lambda wc: BAITS_DICT[wc.bait],
        output:
            accessions=f"{PAIR_DIR}/baits_accessions.txt",
            fasta=f"{PAIR_DIR}/baits.fasta",
        log:
            f"<logs>/make_baits_from_proteinlist/{Path(PAIR_DIR).name}.log",
        wildcard_constraints:
            bait="|".join(Path(b).stem for b in PROTEINLIST_BAITS),
        run:
            copy(input.source, output.accessions)
            accessions = Path(input.source).read_text().splitlines()
            accessions = [x.strip() for x in accessions if x.strip()]
            Path(output.fasta).write_text(fetch_uniprot_fasta(accessions))


if STRUCTURE_BAITS:

    rule make_baits_from_structure:
        input:
            source=lambda wc: BAITS_DICT[wc.bait],
        output:
            accessions=f"{PAIR_DIR}/baits_accessions.txt",
            fasta=f"{PAIR_DIR}/baits.fasta",
        log:
            f"<logs>/make_baits_from_structure/{Path(PAIR_DIR).name}.log",
        wildcard_constraints:
            bait="|".join(Path(b).stem for b in STRUCTURE_BAITS),
        conda:
            f"{ENV_DIR}/base.yaml"
        run:
            uni_ids = get_uniprot_ids_from_struc(input.source)
            fasta = fetch_uniprot_fasta(uni_ids)
            Path(output.accessions).write_text("\n".join(uni_ids) + "\n")
            Path(output.fasta).write_text(fasta)


rule make_preys:
    output:
        preys_fasta=f"{PAIR_DIR}/prey.fasta",
    log:
        f"<logs>/make_preys/{Path(PAIR_DIR).name}.log",
    run:
        with open(output.preys_fasta, "w") as f:
            f.write(fetch_uniprot_fasta([wildcards.prey]))


checkpoint make_manifest:
    input:
        prey_fasta=rules.make_preys.output.preys_fasta,
        baits_fasta=f"{PAIR_DIR}/baits.fasta",
    output:
        manifest=f"{PAIR_DIR}/prediction_manifest.txt",
    params:
        max_len=config.get("total_max_res", 2549),
        database=DATABASE,
    run:
        prey_name, prey_seq = next(iter(fasta_to_dict(input.prey_fasta).items()))
        baits = fasta_to_dict(input.baits_fasta)
        predictions = []
        for bait_name, bait_seq in baits.items():
            name = "__".join(sorted([prey_name, bait_name]))
            if (
                (total_len := len(prey_seq) + len(bait_seq)) > params.max_len
            ) and not Path(f"{params.database}/{name}.tar.zst").exists():
                logger.warning(
                    f"Skipped pair {name} (doesn't exist already and it's too long - {total_len} residues total, max is {params.max_len})"
                )
                continue
            predictions.append(name)
        Path(output.manifest).write_text("\n".join(predictions) + "\n")


rule make_missing_jsons:
    input:
        all_baits=expand(f"{PAIR_DIR}/baits.fasta", bait=BAITS_NAMED, prey=PREYS),
        all_preys=expand(f"{PAIR_DIR}/prey.fasta", bait=BAITS_NAMED, prey=PREYS),
        manifests=expand(
            f"{PAIR_DIR}/prediction_manifest.txt", bait=BAITS_NAMED, prey=PREYS
        ),
    output:
        json_temp=temp(directory(f"<scratch>/missing_jsons")),
    log:
        f"<logs>/make_missing_jsons.log",
    run:
        requested_predictions = []
        seen = set()
        for m in input.manifests:
            pair = Path(m).parent
            for pred in Path(m).read_text().strip().splitlines():
                key = (pair, pred)
                if key not in seen:
                    seen.add(key)
                    requested_predictions.append({"pair": pair, "pred": pred})
        prots = {}
        for p in input.all_preys:
            prots.update(fasta_to_dict(p))
        for b in input.all_baits:
            prots.update(fasta_to_dict(b))
        Path(output.json_temp).mkdir(exist_ok=True)
        for req_pred in requested_predictions:
            pair, name = req_pred["pair"], req_pred["pred"]
            if (pair / f"alphafast_output/{name}/{name}_confidences.json").exists():
                logger.info(f"Skipped pair {pred} (output already exists)")
            else:
                p1, p2 = name.split("__")
                dimer = dimer_template.copy()
                dimer["name"] = name
                dimer["sequences"][0]["protein"]["sequence"] = prots[p1]
                dimer["sequences"][1]["protein"]["sequence"] = prots[p2]
                path = Path(output.json_temp) / f"{name}.json"
                path.write_text(json.dumps(dimer, indent=4))


rule run_all_alphafast:
    input:
        json_temp=rules.make_missing_jsons.output.json_temp,
    output:
        predictions=temp(directory(f"<scratch>/predictions")),
    params:
        modules=modules_for("run_alphafast"),
        af_scriptdir=config["alphafast"]["script_dir"],
        af_container=config["alphafast"]["container"],
        af_backend=config["alphafast"].get("backend", "singularity"),
        af_db_dir=config["alphafast"]["db_dir"],
        af_weights_dir=config["alphafast"]["weights_dir"],
        af_temp_dir=config["alphafast"]["temp_dir"],
        af_jax_dir=config["alphafast"]["jax_compilation_cache_dir"],
        af_gpu_devices=config["alphafast"].get("gpu_devices", "0"),
    shell:
        """
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        mkdir -p {output.predictions}

        export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
        export APPTAINERENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

        {params.modules}

        {params.af_scriptdir}/run_alphafast.sh \
            --input_dir {input.json_temp} \
            --output_dir {output.predictions} \
            --container {params.af_container} \
            --backend {params.af_backend} \
            --db_dir {params.af_db_dir} \
            --weights_dir {params.af_weights_dir} \
            --temp_dir {params.af_temp_dir} \
            --jax_compilation_cache_dir {params.af_jax_dir} \
            --gpu_devices {params.af_gpu_devices}
        """


rule get_predictions:
    input:
        pred=branch(
            lambda wc: wc.prediction in DB_PREDS,
            before_update("<database>/{prediction}.tar.zst"),
            rules.run_all_alphafast.output.predictions,
        ),
    output:
        confidences=temp(
            "<scratch>/prediction_data/{prediction}/{prediction}_confidences.json"
        ),
        summary_confidences=temp(
            "<scratch>/prediction_data/{prediction}/{prediction}_summary_confidences.json"
        ),
        model=temp("<scratch>/prediction_data/{prediction}/{prediction}_model.cif"),
    log:
        "<logs>/get_predictions/{prediction}.log",
    params:
        out_dir=f"{SCRATCH}/prediction_data/{{prediction}}",
    shell:
        r"""
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        if [[ {input.pred:q} == *.tar.zst ]]; then
            tar --extract --zstd \
                --file {input.pred:q} \
                --directory {params.out_dir:q} \
                ./{wildcards.prediction}_confidences.json \
                ./{wildcards.prediction}_model.cif \
                ./{wildcards.prediction}_summary_confidences.json
        else
            ln -s {input.pred:q}/{wildcards.prediction}/{wildcards.prediction}_confidences.json {params.out_dir:q}/{wildcards.prediction}_confidences.json \
                && ln -s {input.pred:q}/{wildcards.prediction}/{wildcards.prediction}_model.cif {params.out_dir:q}/{wildcards.prediction}_model.cif \
                && ln -s {input.pred:q}/{wildcards.prediction}/{wildcards.prediction}_summary_confidences.json {params.out_dir:q}/{wildcards.prediction}_summary_confidences.json
        fi
        """


rule archive_predictions:
    input:
        pred=branch(
            lambda wc: wc.prediction in DB_PREDS,
            [],
            rules.run_all_alphafast.output.predictions,
        ),
    output:
        archive=update("<database>/{prediction}.tar.zst"),
    shell:
        r"""
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        echo "Archiving {input.pred:q}/{wildcards.prediction} as {output.archive:q}"

        if [[ -e {output.archive:q} ]]; then
            echo "Archiving requested for {output.archive:q}, skipping as archive already exists"
            exit 0
        fi

        zst_CLEVEL=19 \
            zst_NBTHREADS={threads} \
            tar --create --zstd \
            --file {output.archive:q} \
            --directory {input.pred:q}/{wildcards.prediction} \
            .
        """


rule map_prediction_to_structure:
    input:
        confs=confidences_for_pair,
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
    output:
        attr_file=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
    log:
        f"<logs>/map_prediction_to_structure/{Path(PAIR_DIR).name}.log",
    conda:
        f"{ENV_DIR}/base.yaml"
    shell:
        r"""
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        python {SCRIPT_DIR}/project_to_structure.py \
            --af3_list {input.confs} \
            --out_file {output.attr_file} \
            --cif_file {input.cif_structure} \
            --n_cpus {threads} \
            --reference_chain {wildcards.prey}
        """


rule get_and_plot_metrics:
    input:
        confs=confidences_for_pair,
        models=models_for_pair,
        summaries=summaries_for_pair,
    output:
        metrics_table=f"{PAIR_DIR}/processed_output/screen_metrics.tsv",
    log:
        f"<logs>/get_and_plot_metrics/{Path(PAIR_DIR).name}.log",
    group:
        "stats"
    conda:
        f"{ENV_DIR}/base.yaml"
    params:
        predictions_dir=f"{SCRATCH}/prediction_data",
        ipsae_dir=f"{PAIR_DIR}/ipSAE",
        out_dir=f"{PAIR_DIR}/processed_output",
        pair=str(Path(PAIR_DIR).name),
        n_files=lambda wc, input: len(input.models),
    shell:
        r"""
        exec > >(tee {log}) 2>&1
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
                python "$ipsae_script" "$jpair" "$model" 10 10 {params.ipsae_dir:q}
            ' _ "{SCRIPT_DIR}/ipsae.py"

        python {SCRIPT_DIR}/ipsae_viz_af3.py \
            --ipsae_dir {params.ipsae_dir:q} \
            --out_dir {params.out_dir:q}

        rm -rf {params.ipsae_dir:q}
        """


rule apply_defattr:
    input:
        cif=lambda wc: BAITS_DICT[wc.bait],
        defattr=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
    output:
        cxs=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.cxs",
    log:
        f"<logs>/apply_defattr/{Path(PAIR_DIR).name}.log",
    params:
        modules=modules_for("apply_defattr"),
        attribute_name=lambda wc: f"{wc.bait}_{wc.prey}",
    shell:
        r"""
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        {params.modules}

        CHIMERAX=$(command -v chimerax || command -v ChimeraX)

        "$CHIMERAX" --nogui <<'EOF'
open {input.cif}
open {input.defattr}
color byattribute r:iCP_{params.attribute_name} #!1 target csab palette 0,#0e8bff:0.5,white:1,red noValueColor #999999
save {output.cxs}
exit
EOF
        """


rule project_screen:
    input:
        confs=confidences_for_pair,
    output:
        out_table=f"{PAIR_DIR}/processed_output/projected_icp.tsv",
    log:
        f"<logs>/project_screen/{Path(PAIR_DIR).name}.log",
    group:
        "stats"
    conda:
        f"{ENV_DIR}/base.yaml"
    params:
        output_dir=f"{PAIR_DIR}/processed_output",
        prey=lambda wc: wc.prey,
    shell:
        r"""
        exec > >(tee {log}) 2>&1
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
            f"{PAIR_DIR}/processed_output/screen_metrics.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
    output:
        summary_table=update(f"{WORK_DIR}/summary_metrics.tsv"),
    log:
        f"<logs>/summarize_screens.log",
    conda:
        f"{ENV_DIR}/base.yaml"
    run:
        logger.info(f"Found {len(input.tables)} tables to concatenate.")
        tables = [pl.scan_csv(t, separator="\t") for t in input.tables]
        summary = pl.concat(tables).sort("ipSAE", descending=True)
        summary.sink_csv(output.summary_table, separator="\t")
        logger.info("Done.")

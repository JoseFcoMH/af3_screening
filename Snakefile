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

BAITS = STRUCTURE_BAITS + PROTEINLIST_BAITS
BAITS_DICT = {Path(b).stem: b for b in BAITS}
BAITS_NAMED = list(BAITS_DICT.keys())

PREYS = Path(config.get("proteins_screened")).read_text().strip().splitlines()

WORK_DIR = config.get("work_dir")

PAIR_DIR = f"{WORK_DIR}/{{bait}}__{{prey}}"
BATCH_DIR = f"{WORK_DIR}/.alphafast_batch"
LOG_DIR = f"{PAIR_DIR}/logs"

logger.info(f"Got PREYS: {PREYS}")
logger.info(f"Got BAITS: {BAITS_NAMED}")


def modules_for(rule_name):
    mods = config.get("modules", {}).get(rule_name, [])
    return "" if not mods else "module load " + " ".join(mods)


def prediction_files_for_pair(bait, prey):
    ckpt = checkpoints.make_manifest.get(bait=bait, prey=prey)
    predictions = Path(ckpt.output.manifest).read_text().strip().splitlines()

    return expand(
        f"{WORK_DIR}/{bait}__{prey}/alphafast_output/{{prediction}}/{{prediction}}_confidences.json",
        prediction=predictions,
    )


def requested_predictions(wc):
    return prediction_files_for_pair(wc.bait, wc.prey)


def all_requested_predictions(wc):
    return [
        f
        for bait in BAITS_NAMED
        for prey in PREYS
        for f in prediction_files_for_pair(bait, prey)
    ]


def jsons_to_predict():
    ckpt = checkpoints.make_missing_jsons.get()
    preds = glob_wildcards(f"{ckpt.output.json_temp}/{{prediction}}.json").prediction
    return expand(f"{ckpt.output.json_temp}/{{prediction}}.json", prediction=preds)


rule all:
    input:
        all_requested_predictions,
        expand(
            f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
            bait=BAITS_NAMED,
            prey=PREYS,
        )
        if structs is not None
        else [],
        expand(
            f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.cxs",
            bait=BAITS_NAMED,
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
            f"{PAIR_DIR}/processed_output/screen_metrics.html",
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
            f"{LOG_DIR}/make_baits_from_proteinlist.log",
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
            f"{LOG_DIR}/make_baits_from_structure.log",
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
        f"{LOG_DIR}/make_preys.log",
    run:
        with open(output.preys_fasta, "w") as f:
            f.write(fetch_uniprot_fasta([wildcards.prey]))


checkpoint make_manifest:
    input:
        prey_fasta=rules.make_preys.output.preys_fasta,
        baits_fasta=f"{PAIR_DIR}/baits.fasta",
    output:
        manifest=temp(f"{PAIR_DIR}/prediction_manifest.txt"),
    params:
        max_len=config.get("total_max_res", 2549),
    run:
        prey_name, prey_seq = next(iter(fasta_to_dict(input.prey_fasta).items()))
        baits = fasta_to_dict(input.baits_fasta)
        predictions = []
        for bait_name, bait_seq in baits.items():
            name = f"{prey_name}__{bait_name}"
            if (total_len := len(prey_seq) + len(bait_seq)) > params.max_len:
                logger.warning(
                    f"Skipped pair {name} (too long - {total_len} residues total, max is {params.max_len})"
                )
                continue
            predictions.append(name)
        Path(output.manifest).write_text("\n".join(predictions) + "\n")


checkpoint make_missing_jsons:
    input:
        all_baits=expand(f"{PAIR_DIR}/baits.fasta", bait=BAITS_NAMED, prey=PREYS),
        all_preys=expand(f"{PAIR_DIR}/prey.fasta", bait=BAITS_NAMED, prey=PREYS),
        manifests=expand(
            f"{PAIR_DIR}/prediction_manifest.txt", bait=BAITS_NAMED, prey=PREYS
        ),
    output:
        json_temp=temp(directory(f"{WORK_DIR}/.affast_input")),
    log:
        f"{WORK_DIR}/make_missing_jsons.log",
    run:
        requested_predictions = []
        for m in input.manifests:
            for pred in Path(m).read_text().strip().splitlines():
                if pred not in requested_predictions:
                    requested_predictions.append(
                        {"pair": Path(m).parent, "pred": pred}
                    )
        preys = {}
        baits = {}
        for p in input.all_preys:
            preys.update(fasta_to_dict(p))
        for b in input.all_baits:
            baits.update(fasta_to_dict(b))
        Path(output.json_temp).mkdir(exist_ok=True)
        for pair, name in requested_predictions.items():
            if (pair / f"alphafast_output/{name}/{name}_confidences.json").exists():
                logger.info(f"Skipped pair {pred} (output already exists)")
            else:
                prey, bait = name.split("__")
                dimer = dimer_template.copy()
                dimer["name"] = name
                dimer["sequences"][0]["protein"]["sequence"] = preys[prey]
                dimer["sequences"][1]["protein"]["sequence"] = baits[bait]
                path = Path(output.json_temp) / f"{name}.json"
                path.write_text(json.dumps(dimer, indent=4))


rule run_all_alphafast:
    input:
        jsons_to_predict,
    output:
        predictions=temp(f"{WORK_DIR}/.predictions_tmp"),
    params:
        json_temp=f"{WORK_DIR}/.affast_input",
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
            --input_dir {params.json_temp} \
            --output_dir {output.predictions} \
            --container {params.af_container} \
            --backend {params.af_backend} \
            --db_dir {params.af_db_dir} \
            --weights_dir {params.af_weights_dir} \
            --temp_dir {params.af_temp_dir} \
            --jax_compilation_cache_dir {params.af_jax_dir} \
            --gpu_devices {params.af_gpu_devices}
        """


rule scatter_predictions:
    input:
        temp_dir=rules.run_all_alphafast.output.predictions,
    output:
        done=f"{PAIR_DIR}/alphafast_output/{{prediction}}/{{prediction}}_confidences.json",
    params:
        predictions_dir=f"{PAIR_DIR}/alphafast_output/{{prediction}}",
    run:
        src_dir = Path(input.temp_dir) / wildcards.prediction
        dst_dir = Path(params.predictions_dir)
        dst_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src_dir, dst_dir)


rule map_prediction_to_structure:
    input:
        requested_predictions,
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
    output:
        attr_file=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
    log:
        f"{LOG_DIR}/map_prediction_to_structure.log",
    conda:
        f"{ENV_DIR}/base.yaml"
    params:
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
    shell:
        """
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        python {SCRIPT_DIR}/AF3screenToStructure.py \
            --af3_dir {params.predictions_dir} \
            --out_file {output.attr_file} \
            --cif_file {input.cif_structure}
        """


rule get_and_plot_metrics:
    input:
        requested_predictions,
    output:
        metrics_table=f"{PAIR_DIR}/processed_output/screen_metrics.tsv",
        metrics_plot=f"{PAIR_DIR}/processed_output/screen_metrics.html",
    log:
        f"{LOG_DIR}/get_and_plot_metrics.log",
    group:
        "stats"
    conda:
        f"{ENV_DIR}/base.yaml"
    params:
        ipsae_dir=f"{PAIR_DIR}/ipSAE",
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
        out_dir=f"{PAIR_DIR}/processed_output",
    shell:
        """
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        mkdir -p {params.ipsae_dir}

        for pair in {params.predictions_dir}/*; do
            echo "Calculating metrics for ${{pair}}"
            pair_name=$(basename $pair)
            jpair=${{pair}}/${{pair_name}}_confidences.json
            model=${{pair}}/${{pair_name}}_model.cif

            if [ -f "$model" ]; then
                python {SCRIPT_DIR}/ipsae.py $jpair $model 10 10
            fi
        done

        mv {params.predictions_dir}/*/*10_10* {params.ipsae_dir}/

        python {SCRIPT_DIR}/ipsae_viz_af3.py \
            --ipsae_dir {params.ipsae_dir} \
            --out_dir {params.out_dir}

        rm -rf {params.ipsae_dir}
        """


rule apply_defattr:
    input:
        cif=lambda wc: BAITS_DICT[wc.bait],
        defattr=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
    output:
        cxs=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.cxs",
    log:
        f"{LOG_DIR}/apply_defattr.log",
    params:
        modules=modules_for("apply_defattr"),
        attribute_name=lambda wc: f"{wc.bait}_{wc.prey}",
    shell:
        """
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
        requested_predictions,
    output:
        out_table=f"{PAIR_DIR}/processed_output/projected_icp.tsv",
    log:
        f"{LOG_DIR}/project_screen.log",
    group:
        "stats"
    conda:
        f"{ENV_DIR}/base.yaml"
    params:
        output_dir=f"{PAIR_DIR}/processed_output",
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
    shell:
        """
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        python {SCRIPT_DIR}/project_screen.py \
            --prediction_dir {params.predictions_dir} \
            --output_dir {params.output_dir} \
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
        summary_table=f"{WORK_DIR}/summary_metrics.tsv",
    log:
        f"{WORK_DIR}/summarize_screens.log",
    conda:
        f"{ENV_DIR}/base.yaml"
    run:
        logger.info(f"Found {len(input.tables)} tables to concatenate.")
        tables = [pl.scan_csv(t, separator="\t") for t in input.tables]
        summary = pl.concat(tables).sort("ipSAE")
        summary.sink_csv(output.summary_table, separator="\t")
        logger.info("Done.")

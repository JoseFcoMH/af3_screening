from snakemake.logging import logger
import logging
from pathlib import Path
import requests
import json
import shutil
from itertools import batched
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

UNIPROT_LOCAL = config.get("uniprot_path", None)
PREDICTION_BATCH_SIZE = int(config.get("prediction_batch_size", 4000))

DB_PATHS = list(Path(config["predictions_database"]).glob("*.tar.zst"))
DB_PREDS = [p.stem.split(".")[0] for p in DB_PATHS]

BAITS = STRUCTURE_BAITS + PROTEINLIST_BAITS
BAITS_DICT = {Path(b).stem: b for b in BAITS}
BAITS_NAMED = list(BAITS_DICT.keys())

PREYS = Path(config.get("proteins_screened")).read_text().strip().splitlines()


logger.info(f"Got PREYS: {PREYS}")
logger.info(f"Got BAITS: {BAITS_NAMED}")
logger.info(f"Baits from protein list: {PROTEINLIST_BAITS}")
logger.info(f"Baits from structure: {STRUCTURE_BAITS}")
logger.info(
    f"Archived predictions in database (showing first 5, out of {len(DB_PATHS)}): {[str(p) for p in DB_PATHS[:5]]}"
)


SCRATCH = config.get("scratch_dir", config["work_dir"] + "/scratch")
DATABASE = config["predictions_database"]
PAIR_NAME = "{bait}__{prey}"
PAIR_DIR = f"<work>/{PAIR_NAME}"
PAIR_OUTPUT = f"{PAIR_DIR}/processed_output"


pathvars:
    database=DATABASE,
    scratch=SCRATCH,
    work=config.get("work_dir"),
    logs="<work>/logs",
    pair_name="{bait}__{prey}",
    pair_dir="<work>/<pair_name>",
    pair_output="<pair_dir>/processed_output",


def modules_for(rule_name):
    mods = config.get("modules", {}).get(rule_name, [])
    return "" if not mods else "module load " + " ".join(mods)


def batch_per_prediction(wc):
    manifest_path = checkpoints.make_batch_manifest.get().output.batch_manifest
    batch_i = (
        pl.scan_csv(manifest_path, separator="\t")
        .filter(pl.col("prediction").eq(wc.prediction))
        .collect()
    )["batch"].to_list()[0]
    return f"<scratch>/predictions/batch_{batch_i}"


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
            f"{PAIR_OUTPUT}/{{bait}}__{{prey}}.defattr",
            bait=[b for b in BAITS_NAMED if BAITS_DICT[b] in STRUCTURE_BAITS],
            prey=PREYS,
        )
        if structs is not None
        else [],
        expand(
            f"{PAIR_OUTPUT}/{{bait}}__{{prey}}.cxs",
            bait=[b for b in BAITS_NAMED if BAITS_DICT[b] in STRUCTURE_BAITS],
            prey=PREYS,
        )
        if structs is not None
        else [],
        expand(
            f"{PAIR_OUTPUT}/screen_metrics.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
        expand(
            f"{PAIR_OUTPUT}/projected_icp.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
        "<work>/summary_metrics.tsv",


rule make_baits_from_proteinlist:
    input:
        source=lambda wc: BAITS_DICT[wc.bait],
    output:
        accessions=f"{PAIR_DIR}/baits_accessions.txt",
        fasta=f"{PAIR_DIR}/baits.fasta",
    log:
        "<logs>/make_baits_from_proteinlist/<pair_name>.log",
    wildcard_constraints:
        bait="|".join(Path(b).stem for b in PROTEINLIST_BAITS),
    run:
        copy(input.source, output.accessions)
        accessions = Path(input.source).read_text().splitlines()
        accessions = [x.strip() for x in accessions if x.strip()]
        Path(output.fasta).write_text(
            fetch_uniprot_fasta(accessions, UNIPROT_LOCAL)
        )


rule make_baits_from_structure:
    input:
        source=lambda wc: BAITS_DICT[wc.bait],
    output:
        accessions=f"{PAIR_DIR}/baits_accessions.txt",
        fasta=f"{PAIR_DIR}/baits.fasta",
    log:
        "<logs>/make_baits_from_structure/<pair_name>.log",
    wildcard_constraints:
        bait="|".join(Path(b).stem for b in STRUCTURE_BAITS),
    conda:
        f"{ENV_DIR}/base.yaml"
    run:
        uni_ids = get_uniprot_ids_from_struc(input.source)
        fasta = fetch_uniprot_fasta(uni_ids, UNIPROT_LOCAL)
        Path(output.accessions).write_text("\n".join(uni_ids) + "\n")
        Path(output.fasta).write_text(fasta)


rule make_preys:
    output:
        preys_fasta="<pair_dir>/prey.fasta",
    log:
        "<logs>/make_preys/<pair_name>.log",
    run:
        with open(output.preys_fasta, "w") as f:
            f.write(fetch_uniprot_fasta([wildcards.prey], UNIPROT_LOCAL))


checkpoint make_manifest:
    input:
        prey_fasta=rules.make_preys.output.preys_fasta,
        baits_fasta="<pair_dir>/baits.fasta",
    output:
        manifest="<pair_dir>/prediction_manifest.txt",
    priority: 10
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


checkpoint make_batch_manifest:
    input:
        all_baits=expand(f"{PAIR_DIR}/baits.fasta", bait=BAITS_NAMED, prey=PREYS),
        all_preys=expand(f"{PAIR_DIR}/prey.fasta", bait=BAITS_NAMED, prey=PREYS),
        manifests=expand(
            f"{PAIR_DIR}/prediction_manifest.txt", bait=BAITS_NAMED, prey=PREYS
        ),
    output:
        batch_manifest=temp("<scratch>/missing_jsons/batched_manifest.tsv"),
        chains=temp("<scratch>/missing_jsons/chains.json"),
    log:
        "<logs>/make_batch_manifest.log",
    params:
        batch_size=PREDICTION_BATCH_SIZE,
    run:
        prots = {}
        for p in input.all_preys + input.all_baits:
            prots.update(fasta_to_dict(p))
        Path(output.chains).write_text(json.dumps(prots, indent=4))
        requested_predictions = set()
        for m in input.manifests:
            for pred in Path(m).read_text().strip().splitlines():
                if spred := pred.strip():
                    requested_predictions.add(spred)
        sorted_predictions = sorted(requested_predictions)
        batches = []
        for batch_i, batch in enumerate(
            batched(sorted_predictions, params.batch_size)
        ):
            batches.append(
                pl.DataFrame({"prediction": batch}).with_columns(
                    pl.lit(batch_i).alias("batch")
                )
            )
        logger.info(
            f"Made {len(batches)} batches of {params.batch_size} predictions (last batch with {len(batches[-1])} predictions)"
        )
        pl.concat(batches).write_csv(output.batch_manifest, separator="\t")


rule make_json_batch:
    input:
        manifest=lambda wc: checkpoints.make_batch_manifest.get().output.batch_manifest,
        chains=lambda wc: checkpoints.make_batch_manifest.get().output.chains,
    output:
        json_temp=temp(
            directory("<scratch>/missing_jsons/batch_{batch}"), group_jobs=True
        ),
    log:
        "<logs>/make_missing_jsons/batch_{batch}.log",
    run:
        json_batch_dir = Path(output.json_temp)
        json_batch_dir.mkdir(exist_ok=True, parents=True)
        with open(input.chains) as f:
            prots = json.load(f)
        manifest = pl.scan_csv(input.manifest, separator="\t")
        prediction_batch = (
            manifest.filter(
                pl.col("batch").cast(pl.Int64).eq(int(wildcards.batch))
            ).collect()
        )["prediction"].to_list()
        for req_pred in prediction_batch:
            p1, p2 = req_pred.split("__")
            dimer = dimer_template.copy()
            dimer["name"] = req_pred
            dimer["sequences"][0]["protein"]["sequence"] = prots[p1]
            dimer["sequences"][1]["protein"]["sequence"] = prots[p2]
            path = json_batch_dir / f"{req_pred}.json"
            path.write_text(json.dumps(dimer, indent=4))


rule run_alphafast_batch:
    input:
        json_temp=rules.make_json_batch.output.json_temp,
    output:
        predictions=temp(
            directory("<scratch>/predictions/batch_{batch}"), group_jobs=True
        ),
    log:
        "<logs>/run_alphafast_batch/batch_{batch}.log",
    priority: 50
    resources:
        af_batch_slots=1,
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


rule archive_predictions:
    input:
        batch_dir=branch(
            lambda wc: wc.prediction in DB_PREDS,
            [],
            batch_per_prediction,
        ),
    output:
        archive=update("<database>/{prediction}.tar.zst"),
    priority: 5
    shell:
        r"""
        exec > >(tee {log}) 2>&1
        set -euo pipefail

        echo "Archiving {input.batch_dir:q}/{wildcards.prediction} as {output.archive:q}"

        if [[ -e {output.archive:q} ]]; then
            echo "Archiving requested for {output.archive:q}, skipping as archive already exists. This should not have happened"
            exit 0
        fi

        zst_CLEVEL=19 \
            zst_NBTHREADS={threads} \
            tar --create --zstd \
            --file {output.archive:q} \
            --directory {input.batch_dir:q}/{wildcards.prediction} \
            .
        """


rule get_predictions:
    input:
        pred=branch(
            lambda wc: wc.prediction in DB_PREDS,
            before_update("<database>/{prediction}.tar.zst"),
            "<database>/{prediction}.tar.zst",
        ),
    output:
        confidences=temp(
            "<scratch>/prediction_data/{prediction}/{prediction}_confidences.json",
            group_jobs=True,
        ),
        summary_confidences=temp(
            "<scratch>/prediction_data/{prediction}/{prediction}_summary_confidences.json",
            group_jobs=True,
        ),
        model=temp(
            "<scratch>/prediction_data/{prediction}/{prediction}_model.cif",
            group_jobs=True,
        ),
    log:
        "<logs>/get_predictions/{prediction}.log",
    priority: -5
    params:
        out_dir=lambda wc, output: str(Path(output.confidences).parent),
    shell:
        r"""
        exec > >(tee {log}) 2>&1
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
        f"{ENV_DIR}/base.yaml"
    params:
        ipsae_dir=lambda wc: f"{SCRATCH}/{wc.bait}__{wc.prey}/ipSAE",
        out_dir=lambda wc, output: str(Path(output.metrics_table).parent),
        pair=lambda wc: f"{wc.bait}__{wc.prey}",
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


rule map_prediction_to_structure:
    input:
        confs=confidences_for_pair,
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
    output:
        attr_file="<pair_output>/{bait}__{prey}.defattr",
    log:
        "<logs>/map_prediction_to_structure/<pair_name>.log",
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


rule apply_defattr:
    input:
        cif=lambda wc: BAITS_DICT[wc.bait],
        defattr="<pair_output>/{bait}__{prey}.defattr",
    output:
        cxs="<pair_output>/{bait}__{prey}.cxs",
    log:
        "<logs>/apply_defattr/<pair_name>.log",
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
        out_table="<pair_output>/projected_icp.tsv",
    log:
        "<logs>/project_screen/<pair_name>.log",
    group:
        "stats"
    conda:
        f"{ENV_DIR}/base.yaml"
    params:
        output_dir=lambda wc, output: str(Path(output.out_table).parent),
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
            f"{PAIR_OUTPUT}/screen_metrics.tsv",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
    output:
        summary_table=update("<work>/summary_metrics.tsv"),
    log:
        "<logs>/summarize_screens.log",
    conda:
        f"{ENV_DIR}/base.yaml"
    run:
        logger.info(f"Found {len(input.tables)} tables to concatenate.")
        tables = [pl.scan_csv(t, separator="\t") for t in input.tables]
        summary = pl.concat(tables).sort("ipSAE", descending=True).unique("pair")
        summary.sink_csv(output.summary_table, separator="\t")
        logger.info("Done.")

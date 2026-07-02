from snakemake.logging import logger
import logging
from pathlib import Path
import requests
import json
from shutil import copy, copytree

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

LOG_DIR = f"{PAIR_DIR}/logs"

logger.info(f"Got PREYS: {PREYS}")
logger.info(f"Got BAITS: {BAITS_NAMED}")


def modules_for(rule_name):
    mods = config.get("modules", {}).get(rule_name, [])
    return "" if not mods else "module load " + " ".join(mods)


rule all:
    input:
        expand(
            f"{PAIR_DIR}/alphafast_output/.done",
            bait=BAITS_NAMED,
            prey=PREYS,
        )
        if False
        else [],
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


rule make_input_json:
    input:
        all_baits=expand(f"{PAIR_DIR}/baits.fasta", bait=BAITS_NAMED, prey=PREYS),
        all_preys=expand(f"{PAIR_DIR}/prey.fasta", bait=BAITS_NAMED, prey=PREYS),
    output:
        json_temp=temp(directory(f"{WORK_DIR}/.affast_input")),
    log:
        f"{WORK_DIR}/make_input_json.log",
    params:
        max_len=config.get("total_max_res", 2549),
    run:
        preys = {}
        baits = {}
        for p in input.all_preys:
            preys.update(fasta_to_dict(p))
        for b in input.all_baits:
            baits.update(fasta_to_dict(b))
        Path(output.json_temp).mkdir(exist_ok=True)
        for prey_name, prey_seq in preys.items():
            for bait_name, bait_seq in baits.items():
                name = f"{prey_name}__{bait_name}"
                if (total_len := len(prey_seq) + len(bait_seq)) > params.max_len:
                    logger.warning(
                        f"Skipped pair {name} (too long - {total_len} residues total, max is {params.max_len})"
                    )
                    continue
                if (
                    len(
                        list(
                            Path(WORK_DIR).glob(f"*/alphafast_output/{name}/*.cif")
                        )
                    )
                    > 0
                ):
                    logger.info(
                        f"Skipped pair {name} (output structure already exists)"
                    )
                    continue
                dimer = dimer_template.copy()
                dimer["name"] = name
                dimer["sequences"][0]["protein"]["sequence"] = prey_seq
                dimer["sequences"][1]["protein"]["sequence"] = bait_seq
                path = Path(output.json_temp) / f"{name}.json"
                path.write_text(json.dumps(dimer, indent=4))


rule run_all_alphafast:
    input:
        json_temp=rules.make_input_json.output.json_temp,
    output:
        predictions=temp(directory(f"{WORK_DIR}/.predictions_tmp")),
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


rule scatter_predictions:
    input:
        temp_dir=rules.run_all_alphafast.output.predictions,
        bait_accessions=f"{PAIR_DIR}/baits_accessions.txt",
    output:
        done=touch(f"{PAIR_DIR}/alphafast_output/.done"),
    params:
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
    run:
        Path(params.predictions_dir).mkdir(exist_ok=True)
        baits = Path(input.bait_accessions).read_text().strip().splitlines()
        for pred in Path(input.temp_dir).glob("*"):
            pred_names = pred.name.split("_")
            for bait in baits:
                if (bait in pred_names) and (wildcards.prey in pred_names):
                    logger.info(
                        f"Copying {pred.name} into {params.predictions_dir}"
                    )
                    copytree(
                        pred,
                        Path(params.predictions_dir) / pred.name,
                        dirs_exist_ok=True,
                    )


rule map_prediction_to_structure:
    input:
        rules.scatter_predictions.output.done,
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
        rules.scatter_predictions.output.done,
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
        rules.scatter_predictions.output.done,
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

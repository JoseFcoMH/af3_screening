from snakemake.logging import logger
from pathlib import Path
import requests
import json
from shutil import copy, rmtree
import yaml

from load_utils import (
    fasta_to_dict,
    fetch_uniprot_fasta,
    get_uniprot_ids_from_struc,
    dimer_template,
)


SCRIPT_DIR = workflow.source_path("scripts")
ENV_DIR = workflow.source_path("envs")
CONFIG_DIR = workflow.source_path("config")


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

logger.info(f"Got PREYS: {PREYS}")
logger.info(f"Got BAITS: {BAITS_NAMED}")


def modules_for(rule_name):
    mods = config.get("modules", {}).get(rule_name, [])
    return "" if not mods else "module load " + " ".join(mods)


rule all:
    input:
        expand(
            f"{PAIR_DIR}/alphafast_output",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
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
    run:
        with open(output.preys_fasta, "w") as f:
            f.write(fetch_uniprot_fasta([wildcards.prey]))


rule make_input_json:
    input:
        baits_fasta=f"{PAIR_DIR}/baits.fasta",
        preys_fasta=f"{PAIR_DIR}/prey.fasta",
    output:
        json_dir=directory(f"{PAIR_DIR}/data_pipeline_input"),
    params:
        max_len=config.get("total_max_res", 2549),
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
    run:
        preys = fasta_to_dict(input.preys_fasta)
        baits = fasta_to_dict(input.baits_fasta)

        if Path(output.json_dir).exists():
            rmtree(output.json_dir)
        Path(output.json_dir).mkdir(parents=True)

        for prey_name, prey_seq in preys.items():
            for bait_name, bait_seq in baits.items():

                name = f"{prey_name}__{bait_name}"

                if (total_len := len(prey_seq) + len(bait_seq)) > params.max_len:
                    logger.warning(
                        f"Skipped pair {name} (too long - {total_len} residues total, max is {params.max_len})"
                    )
                    continue
                if (Path(params.predictions_dir) / f"{name}/{name}_model.cif").exists():
                    logger.info(
                        f"Skipped pair {name} (output structure already exists)"
                    )
                    continue

                dimer = dimer_template.copy()

                dimer["name"] = name
                dimer["sequences"][0]["protein"]["sequence"] = prey_seq
                dimer["sequences"][1]["protein"]["sequence"] = bait_seq

                path = Path(output.json_dir) / f"{name}.json"
                path.write_text(json.dumps(dimer, indent=4))


rule run_alphafast:
    input:
        json_dir=f"{PAIR_DIR}/data_pipeline_input",
    output:
        predictions_dir=directory(f"{PAIR_DIR}/alphafast_output"),
    params:
        modules=modules_for("run_alphafast"),
    shell:
        """
    if [ -z "$(ls -A "{input.json_dir}")" ]; then
      echo 'No json to process at {input.json_dir}, skipping alignment and inference'
      exit 0
    fi

    export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
    export APPTAINERENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

    {params.modules}

    {config['alphafast']['script_dir']}/run_alphafast.sh \
      --input_dir {input.json_dir} \
      --output_dir {output.predictions_dir} \
      --container {config['alphafast']['container']} \
      --db_dir {config['alphafast']['db_dir']} \
      --weights_dir {config['alphafast']['weights_dir']} \
      --temp_dir {config['alphafast']['temp_dir']} \
      --jax_compilation_cache_dir {config['alphafast']['temp_dir']}
    """


rule map_prediction_to_structure:
    input:
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
    output:
        attr_file=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
    conda:
        f"{ENV_DIR}/base.yaml"
    shell:
        """
      python {SCRIPT_DIR}/AF3screenToStructure.py \
          --af3_dir {input.predictions_dir} \
          --out_file {output.attr_file} \
          --cif_file {input.cif_structure}  \
    """


rule get_and_plot_metrics:
    input:
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
    output:
        metrics_table=f"{PAIR_DIR}/processed_output/screen_metrics.tsv",
        metrics_plot=f"{PAIR_DIR}/processed_output/screen_metrics.html",
    params:
        ipsae_dir=f"{PAIR_DIR}/ipSAE",
        out_dir=f"{PAIR_DIR}/processed_output",
    conda:
        f"{ENV_DIR}/base.yaml"
    group:
        "stats"
    shell:
        """ 
        mkdir -p {params.ipsae_dir} 

        for pair in {input.predictions_dir}/*; do
            echo "Calculating metrics for ${{pair}}"
            pair_name=$(basename $pair)
            jpair=${{pair}}/${{pair_name}}_confidences.json
            model=${{pair}}/${{pair_name}}_model.cif
            python {SCRIPT_DIR}/ipsae.py $jpair $model 10 10
        done

        mv {input.predictions_dir}/*/*10_10* {params.ipsae_dir}/

        python {SCRIPT_DIR}/ipsae_viz_af3.py \
          --ipsae_dir {params.ipsae_dir} \
          --out_dir {params.out_dir}
    """


rule apply_defattr:
    input:
        cif_structure=lambda wc: BAITS_DICT[wc.bait],
        attr_file=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.defattr",
    output:
        attr_structure=f"{PAIR_DIR}/processed_output/{{bait}}__{{prey}}.cxs",
    params:
        modules=modules_for("apply_defattr"),
    shell:
        """
        {params.modules}

        ChimeraX --nogui <<'EOF'
        open {input.cif_structure}
        open {input.attr_file}
        color byattribute r:iCP_{wildcards.bait}_{wildcards.prey} #!1 target csab palette 0,#0e8bff:0.5,white:1,red noValueColor #999999
        save {output.attr_structure}
        exit
        EOF
    """


rule project_screen:
    input:
        predictions_dir=f"{PAIR_DIR}/alphafast_output",
    output:
        out_table=f"{PAIR_DIR}/processed_output/projected_icp.tsv",
    params:
        output_dir=f"{PAIR_DIR}/processed_output",
    group:
        "stats"
    conda:
        f"{ENV_DIR}/base.yaml"
    shell:
        """
        python {SCRIPT_DIR}/project_screen.py \
            --prediction_dir {input.predictions_dir} \
            --output_dir {params.output_dir} \
            --n_cores {threads}
    """

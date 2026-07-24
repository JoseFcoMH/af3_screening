import json
import polars as pl
from pathlib import Path
from scripts.load_utils import dimer_template


rule make_json_batch:
    input:
        manifests=lambda wc: manifests_for_batches(),
        fasta=ancient(rules.make_master_fasta.output.fasta),
    output:
        json_temp=temp(
            directory("<scratch>/missing_jsons/batch_{batch}"), group_jobs=True
        ),
    log:
        "<logs>/make_missing_jsons/batch_{batch}.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            json_batch_dir = Path(output.json_temp)
            json_batch_dir.mkdir(exist_ok=True, parents=True)
            prots = fasta_to_dict(input.fasta)
            manifest = batch_manifests(input.manifests)
            prediction_batch = (
                manifest.filter(
                    pl.col("batch").cast(pl.Int64).eq(int(wildcards.batch))
                )
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
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        mkdir -p {output.predictions:q}

        export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
        export APPTAINERENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

        {params.modules}

        {params.af_scriptdir}/run_alphafast.sh \
            --input_dir {input.json_temp:q} \
            --output_dir {output.predictions:q} \
            --container {params.af_container:q} \
            --backend {params.af_backend:q} \
            --db_dir {params.af_db_dir:q} \
            --weights_dir {params.af_weights_dir:q} \
            --temp_dir {params.af_temp_dir:q} \
            --jax_compilation_cache_dir {params.af_jax_dir:q} \
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
    log:
        "<logs>/archive_predictions/{prediction}.log",
    priority: 20
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
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

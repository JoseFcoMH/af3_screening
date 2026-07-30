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


rule setup_alphafast:
    output:
        af_run_script=str(
            Path(config.get("af_repo_dir", "<scratch>/alphafast/repo"))
            / "scripts/run_alphafast.sh"
        ),
        mmseqs_sentinel=str(
            Path(config.get("af_mmseqs_db", "<scratch>/alphafast/mmseqs"))
            / "mmseqs.ready"
        ),
    log:
        "<logs>/setup_alphafast.log",
    conda:
        ENVS["udocker"]
    params:
        af_repo=config.get("af_repo", "https://github.com/JoseFcoMH/alphafast.git"),
        repo_parent=lambda wc, output: str(Path(output.af_run_script).parent.parent),
        mmseqs_dir=lambda wc, output: str(Path(output.mmseqs_sentinel).parent),
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        echo "=== Copying repo ==="
        cd {params.repo_parent:q}
        git clone {params.af_repo:q}

        echo "=== Setting up mmseqs database ==="
        ./scripts/setup_databases.sh {params.mmseqs_dir:q} --protein-only
        touch {output.mmseqs_sentinel:q}
        """


rule run_alphafast_batch:
    input:
        af_run_script=ancient(rules.setup_alphafast.output.af_run_script),
        mmseqs_ready=ancient(rules.setup_alphafast.output.mmseqs_sentinel),
        json_temp=rules.make_json_batch.output.json_temp,
    output:
        predictions=temp(
            directory("<scratch>/predictions/batch_{batch}"), group_jobs=True
        ),
    log:
        "<logs>/run_alphafast_batch/batch_{batch}.log",
    priority: 50
    conda:
        ENVS["udocker"]
    resources:
        af_batch_slots=1,
        n_gpu=1,
    params:
        weights_dir=config["af_weights_dir"],
        container=config.get("af_container", "romerolabduke/alphafast:latest"),
        backend=config.get("af_backend", "udocker"),
        temp_dir=str(Path(SCRATCH) / "alphafast"),
        jax_dir=str(Path(SCRATCH) / "jax"),
        mmseqs_db=lambda wc, input: str(Path(input.mmseqs_ready).parent.parent),
        gpu_devices=lambda wc, resources: ",".join(
            str(x) for x in range(int(resources.n_gpu))
        ),
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        mkdir -p {output.predictions:q}

        {input.af_run_script:q} \
            --input_dir {input.json_temp:q} \
            --output_dir {output.predictions:q} \
            --container {params.container:q} \
            --backend {params.backend} \
            --db_dir {params.mmseqs_db:q} \
            --weights_dir {params.weights_dir:q} \
            --temp_dir {params.temp_dir:q} \
            --jax_compilation_cache_dir {params.jax_dir:q} \
            --gpu_devices {params.gpu_devices}
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
    priority: 100
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

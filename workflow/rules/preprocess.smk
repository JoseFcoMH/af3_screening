from itertools import batched
from functools import cache
from types import MappingProxyType
from pathlib import Path
import polars as pl

from scripts.load_utils import (
    extract_fasta_records,
    index_fasta,
    fasta_to_dict,
    fetch_uniprot_fasta,
    get_uniprot_ids_from_struc,
)


def batch_manifests(
    manifests: list[Path | str],
    batch_size: int = PREDICTION_BATCH_SIZE,
    done: Path | str = DB_PREDS,
) -> pl.DataFrame:
    requested_predictions = set()
    for m in manifests:
        for pred in Path(m).read_text().strip().splitlines():
            if pred not in done:
                requested_predictions.add(pred)
    sorted_predictions = sorted(requested_predictions)
    batches = []
    for batch_i, batch in enumerate(batched(sorted_predictions, batch_size)):
        batches.append(
            pl.DataFrame({"prediction": batch}).with_columns(
                pl.lit(batch_i).alias("batch")
            )
        )
    logger.info(
        f"Made {len(batches)} batches of {len(batches[0])} predictions (last batch with {len(batches[-1])} predictions)"
    )
    return pl.concat(batches)


def manifests_for_batches() -> list[str]:
    manifest_path = checkpoints.make_project_manifest.get().output.project_manifest
    manifest = pl.scan_csv(manifest_path, separator="\t")
    manifest_missing_preds = (
        manifest.filter(~pl.col("pair").is_in(DB_PREDS))
        .select("from_manifest")
        .unique("from_manifest", maintain_order=True)
    )
    batch_manifests = manifest_missing_preds.collect()["from_manifest"].to_list()
    return batch_manifests


@cache
def make_batch_mapping() -> MappingProxyType:
    manifest_paths = manifests_for_batches()
    batched_manifests = batch_manifests(manifest_paths)
    mapping = {
        row["prediction"]: row["batch"]
        for row in batched_manifests.iter_rows(named=True)
    }
    return MappingProxyType(mapping)


def batch_per_prediction(wc) -> str:
    manifest_mapping = make_batch_mapping()
    batch_i = manifest_mapping[wc.prediction]
    return f"<scratch>/predictions/batch_{batch_i}"


@cache
def read_pair_manifest(bait: str, prey: str) -> tuple[str]:
    manifest_path = checkpoints.make_project_manifest.get().output.project_manifest
    manifest = pl.scan_csv(manifest_path, separator="\t")
    pair_manifest = manifest.filter(pl.col("prey").eq(prey) & pl.col("bait").eq(bait))
    pair_predictions = pair_manifest.collect()["prediction"].to_list()
    return tuple(pair_predictions)


def confidences_for_pair(wc) -> list[str]:
    bait, prey = wc.bait, wc.prey
    predictions = read_pair_manifest(bait=bait, prey=prey)
    return expand(
        "<pair_scratch>/{prediction}/{prediction}_confidences.json",
        prediction=predictions,
    )


def models_for_pair(wc) -> list[str]:
    bait, prey = wc.bait, wc.prey
    predictions = read_pair_manifest(bait=bait, prey=prey)
    return expand(
        "<pair_scratch>/{prediction}/{prediction}_model.cif",
        prediction=predictions,
    )


def summaries_for_pair(wc) -> list[str]:
    bait, prey = wc.bait, wc.prey
    predictions = read_pair_manifest(bait=bait, prey=prey)
    return expand(
        "<pair_scratch>/{prediction}/{prediction}_summary_confidences.json",
        prediction=predictions,
    )


def predictions_for_pair(bait, prey) -> list[str]:
    predictions = read_pair_manifest(bait=bait, prey=prey)
    return expand(
        "<database>/{prediction}.tar.zst",
        prediction=predictions,
    )


def all_predictions(wc) -> list[str]:
    preds = set()
    for bait in BAITS_NAMED:
        for prey in PREYS:
            pred_names = read_pair_manifest(bait=bait, prey=prey)
            preds.update([f"{DATABASE}/{n}.tar.zst" for n in pred_names])
    return sorted(preds)


if PROTEINLIST_BAITS:

    rule make_baits_from_proteinlist:
        input:
            source=lambda wc: BAITS_DICT[wc.bait],
        output:
            accessions=f"{PAIR_DIR}/baits_accessions.txt",
        log:
            "<logs>/make_baits_from_proteinlist/<pair_name>.log",
        wildcard_constraints:
            bait="|".join(Path(b).stem for b in PROTEINLIST_BAITS),
        shell:
            r"""
            exec > >(tee {log:q}) 2>&1
            set -euo pipefail

            sort -u {input.source:q} >{output.accessions:q}
            echo "Number of unique IDs in {output.accessions:q}: $(wc -l {output.accessions:q})"
            """


if STRUCTURE_BAITS:

    rule make_baits_from_structure:
        input:
            source=lambda wc: BAITS_DICT[wc.bait],
        output:
            accessions=f"{PAIR_DIR}/baits_accessions.txt",
        log:
            "<logs>/make_baits_from_structure/<pair_name>.log",
        wildcard_constraints:
            bait="|".join(Path(b).stem for b in STRUCTURE_BAITS),
        run:
            with rule_logger(logfile=log[0], rule_name=rule) as logger:
                uni_ids = get_uniprot_ids_from_struc(input.source)
                Path(output.accessions).write_text("\n".join(uni_ids) + "\n")


rule make_uniprot_index:
    input:
        database=UNIPROT_LOCAL,
    output:
        index=str(Path(UNIPROT_LOCAL).parent / f"{Path(UNIPROT_LOCAL).stem}.idx"),
    log:
        "<logs>/make_uniprot_index.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            index_fasta(fasta_path=input.database, idx_path=output.index)


rule make_master_fasta:
    input:
        bait_lists=expand(
            f"{PAIR_DIR}/baits_accessions.txt", prey=PREYS, bait=BAITS_NAMED
        ),
        uniprot_db=ancient(UNIPROT_LOCAL),
        uniprot_idx=ancient(rules.make_uniprot_index.output.index),
    output:
        accessions="<work>/all_accessions.txt",
        fasta="<work>/all_chains.fasta",
        fasta_index="<work>/all_chains.idx",
    log:
        "<logs>/make_master_fasta.log",
    params:
        prey_accessions=PREYS,
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            all_accessions = set(params.prey_accessions)
            for blist in input.bait_lists:
                baccesssions = Path(blist).read_text().strip().splitlines()
                all_accessions.update(baccesssions)
            sorted_accessions = sorted(all_accessions)
            Path(output.accessions).write_text("\n".join(sorted_accessions))
            fasta = fetch_uniprot_fasta(sorted_accessions, input.uniprot_idx)
            Path(output.fasta).write_text(fasta)
            index_fasta(fasta_path=output.fasta, idx_path=output.fasta_index)


rule make_prey_fasta_per_pair:
    input:
        fasta_index=ancient(rules.make_master_fasta.output.fasta_index),
    output:
        prey_fasta="<pair_dir>/prey.fasta",
    log:
        "<logs>/make_prey_fasta_per_pair/<pair_name>.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            _, prey_fasta = extract_fasta_records(
                index_path=input.fasta_index, accessions=[wildcards.prey]
            )
            Path(output.prey_fasta).write_text(prey_fasta)


rule make_baits_fasta_per_pair:
    input:
        fasta_index=ancient(rules.make_master_fasta.output.fasta_index),
        baits_accessions="<pair_dir>/baits_accessions.txt",
    output:
        baits_fasta="<pair_dir>/baits.fasta",
    log:
        "<logs>/make_baits_fasta_per_pair/<pair_name>.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            baits = Path(input.baits_accessions).read_text().strip().splitlines()
            _, baits_fasta = extract_fasta_records(
                index_path=input.fasta_index, accessions=baits
            )
            Path(output.baits_fasta).write_text(baits_fasta)


rule make_manifest:
    input:
        prey_fasta="<pair_dir>/prey.fasta",
        baits_fasta="<pair_dir>/baits.fasta",
    output:
        manifest="<pair_dir>/prediction_manifest.txt",
    log:
        "<logs>/make_manifests/<pair_name>.log",
    params:
        max_len=config.get("total_max_res", 2549),
        database=DATABASE,
        done_preds=DB_PREDS,
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            prey_name, prey_seq = next(
                iter(fasta_to_dict(input.prey_fasta).items())
            )
            baits = fasta_to_dict(input.baits_fasta)
            predictions = []
            for bait_name, bait_seq in baits.items():
                name = "__".join(sorted([prey_name, bait_name]))
                if (
                    (total_len := len(prey_seq) + len(bait_seq)) > params.max_len
                ) and name not in params.done_preds:
                    logger.warning(
                        f"Skipped pair {name} (doesn't already exist and it's too long - {total_len} residues total, max is {params.max_len})"
                    )
                    continue
                predictions.append(name)
            logger.info(f"Wrote {len(predictions)} entries in manifest.")
            Path(output.manifest).write_text("\n".join(predictions))


checkpoint make_project_manifest:
    input:
        manifests=expand(
            f"{PAIR_DIR}/prediction_manifest.txt",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
    output:
        project_manifest="<work>/project_manifest.tsv",
    log:
        "<logs>/make_project_manifest.log",
    priority: 100
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            logger.info(f"Found {len(input.manifests)} manifests to concatenate.")
            manifest_tables = []
            for mani in input.manifests:
                bait, prey = Path(mani).parent.name.split("__")
                mani_content = Path(mani).read_text().strip().splitlines()
                mani_table = pl.DataFrame(
                    {"prediction": mani_content}
                ).with_columns(
                    pl.lit(str(mani)).alias("from_manifest"),
                    pl.lit(bait).alias("bait"),
                    pl.lit(prey).alias("prey"),
                    pl.lit(f"{bait}__{prey}").alias("pair"),
                )
                manifest_tables.append(mani_table)
            project_manifest = pl.concat(manifest_tables).sort(
                ["from_manifest", "pair", "prediction"]
            )
            project_manifest.write_csv(output.project_manifest, separator="\t")
            logger.info(
                f"Found {len(project_manifest)} requested predictions, {len(project_manifest.unique('prediction'))} of them unique."
            )

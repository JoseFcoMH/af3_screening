from itertools import batched
from functools import cache
from types import MappingProxyType
from pathlib import Path
import polars as pl


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

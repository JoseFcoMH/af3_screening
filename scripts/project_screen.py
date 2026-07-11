#!/usr/bin/env python3

import argparse
from pathlib import Path
import json
import numpy as np
import altair as alt
import polars as pl
import logging
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


@dataclass
class DimerPrediction:
    prediction: Path

    def __post_init__(self):
        with open(self.prediction / f"{self.name}_confidences.json") as f:
            data = json.load(f)
            self.chain_length = {
                "A": sum(t == "A" for t in data["token_chain_ids"]),
                "B": sum(t == "B" for t in data["token_chain_ids"]),
            }
            probs = np.array(data["contact_probs"])
            interface1 = probs[: self.chain_length["A"], self.chain_length["A"] :].max(
                1
            )
            interface2 = probs[self.chain_length["A"] :, : self.chain_length["A"]].max(
                0
            )
            self.icp = np.stack([interface1, interface2]).max(0)

    @property
    def name(self):
        return self.prediction.stem

    @property
    def names(self):
        name_a, name_b = self.name.split("__")
        return {"A": name_a, "B": name_b}

    @property
    def icp_table(self) -> pl.DataFrame:
        len_a = self.chain_length["A"]
        return pl.DataFrame(
            {
                "residue": list(range(1, len_a + 1)),
                "icp": self.icp,
                "name_A": [self.names["A"]] * len_a,
                "name_B": [self.names["B"]] * len_a,
            }
        )


def load_prediction(prediction: Path) -> pl.DataFrame:
    return DimerPrediction(prediction).icp_table


def find_predictions(prediction_dir: Path) -> list[Path]:
    predictions = []
    logger.info(f"Scanning {prediction_dir}...")
    for pred in prediction_dir.glob("*"):
        if len(list(pred.glob("*.cif"))) == 1:
            predictions.append(pred)
        else:
            logger.info(f"Skipped {pred} (no structure found in dir)")
    logger.info(f"Found {len(predictions)} usable predictions")
    return predictions


def make_chart(data: pl.DataFrame) -> alt.FacetChart:
    line = (
        alt.Chart(data)
        .mark_line()
        .encode(
            x=alt.X("residue:Q", title="Residue #"),
            y=alt.Y(
                "icp:Q",
                title="Contact probability",
                scale=alt.Scale(domain=[0, 1]),
            ),
            color=alt.Color("name_B:N", title="Interactor"),
        )
    )

    points = (
        alt.Chart(data)
        .transform_filter("datum.icp > 0.4")
        .mark_point(size=40)
        .encode(
            x="residue:Q",
            y="icp:Q",
            color="name_B:N",
            tooltip=[
                alt.Tooltip("name_B:N", title="Interactor"),
                alt.Tooltip("residue:Q", title="Residue"),
                alt.Tooltip("icp:Q", title="Contact probability", format=".2f"),
            ],
        )
    )

    chart = (
        alt.layer(line, points)
        .properties(
            height=400,
            width=600,
        )
        .interactive()
        .facet(row=alt.Row("name_A:N", title="Prey"))
    )
    return chart


def main():

    parser = argparse.ArgumentParser()

    # parser.add_argument("--prediction_dir", type=Path)
    parser.add_argument("--predictions_list", nargs="+", type=Path)
    parser.add_argument("--output_dir", type=Path)
    parser.add_argument("--n_cores", type=int, default=1)

    args = parser.parse_args()
    # prediction_dir = args.prediction_dir
    output_dir = args.output_dir

    out_table = output_dir / "projected_icp.tsv"
    out_html = output_dir / "projected_icp.html"
    out_png = output_dir / "projected_icp.png"

    # predictions = find_predictions(prediction_dir)
    predictions = args.predictions_list

    logger.info("Started processing predictions.")
    with ProcessPoolExecutor(args.n_cores) as executor:
        results = list(executor.map(load_prediction, predictions))
    results_table = pl.concat(results)
    results_table.write_csv(out_table, separator="\t")

    logger.info(f"Done. Saved results as {out_table}. Started plotting.")

    hit_genes = results_table.filter(pl.col("icp").gt(0.6)).unique(["name_A", "name_B"])
    filtered_data = results_table.join(hit_genes, on=["name_A", "name_B"], how="semi")

    if len(filtered_data) == 0:
        logger.info("No hits, so nothing to plot")
    else:
        logger.info(f"Kept {len(filtered_data.unique('name_B'))} hits to plot.")
        chart = make_chart(filtered_data)
        chart.save(out_html, inline=True, embed_options={"renderer": "canvas"})
        chart.save(out_png, scale_factor=3.0)
        logger.info("Saved plots.")

    logger.info("All done.")


if __name__ == "__main__":
    main()

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
from itertools import repeat

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


@dataclass
class DimerPrediction:
    prediction: Path
    reference_chain: str

    def __post_init__(self):
        with open(self.prediction / f"{self.name}_confidences.json") as f:
            data = json.load(f)
            self.probs = np.array(data["contact_probs"]).astype(np.float32)
            self.chain_length = {
                "A": sum(t == "A" for t in data["token_chain_ids"]),
                "B": sum(t == "B" for t in data["token_chain_ids"]),
            }
            self.project_on = (
                "A" if self.name.split("_")[0] == self.reference_chain else "B"
            )
            logging.debug(
                f"Will project on chain {self.project_on} aka {self.names[self.project_on]}"
            )
            self.icp = self.project_icp(self.project_on)

    @property
    def name(self) -> str:
        return self.prediction.stem

    @property
    def names(self) -> dict[str, str]:
        name_a, name_b = self.name.split("__")
        return {"A": name_a, "B": name_b}

    @property
    def icp_table(self) -> pl.DataFrame:
        len_projection = (
            self.chain_length["A"] if self.project_on == "A" else self.chain_length["B"]
        )

        return pl.DataFrame(
            {
                "residue": list(range(1, len_projection + 1)),
                "icp": self.icp,
                "pair": self.name,
                "name_A": self.names["A"],
                "name_B": self.names["B"],
                "project_on": self.project_on,
            }
        )

    def project_icp(self, on="A") -> np.ndarray:
        interface1 = self.probs[: self.chain_length["A"], self.chain_length["A"] :].max(
            1 if on == "A" else 0
        )
        interface2 = self.probs[self.chain_length["A"] :, : self.chain_length["A"]].max(
            0 if on == "A" else 1
        )
        return np.stack([interface1, interface2]).max(0)


def load_prediction(prediction: Path, project_on: str) -> pl.DataFrame:
    return DimerPrediction(prediction, reference_chain=project_on).icp_table


def make_chart(data: pl.DataFrame) -> alt.FacetChart:
    data = data.with_columns(
        pl.when(pl.col("project_on") == "A")
        .then(pl.col("name_B"))
        .otherwise(pl.col("name_A"))
        .alias("projected"),
        pl.col("icp")
        .max()
        .over("pair")
        .rank(method="dense", descending=True)
        .alias("pair_rank"),
    ).filter(pl.col("pair_rank").lt(13))

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
            color=alt.Color("projected:N", title="Interactor"),
        )
    )

    points = (
        alt.Chart(data)
        .transform_filter("datum.icp > 0.6")
        .mark_point(size=30)
        .encode(
            x="residue:Q",
            y="icp:Q",
            color="projected:N",
            tooltip=[
                alt.Tooltip("projected:N", title="Interactor"),
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
    )
    return chart


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--predictions_list", type=Path)
    parser.add_argument("--output_dir", type=Path)
    parser.add_argument("--prey", type=str)
    parser.add_argument("--n_cores", type=int, default=1)

    args = parser.parse_args()

    output_dir = args.output_dir
    prey = args.prey

    out_table = output_dir / "projected_icp.tsv"
    out_html = output_dir / "projected_icp.html"
    out_png = output_dir / "projected_icp.png"

    predictions = [
        Path(p).parent
        for p in Path(args.predictions_list).read_text().strip().splitlines()
    ]

    logger.info("Started processing predictions.")
    with ProcessPoolExecutor(args.n_cores) as executor:
        results = list(executor.map(load_prediction, predictions, repeat(prey)))
    results_table = pl.concat(results)
    results_table.write_csv(out_table, separator="\t")

    logger.info(f"Done. Saved results as {out_table}. Started plotting.")

    hit_genes = results_table.filter(pl.col("icp").gt(0.4)).unique(["name_A", "name_B"])
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

#!/usr/bin/env python3

import argparse
from pathlib import Path
import json
import numpy as np
import polars as pl
import altair as alt
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


def make_chart(
    data: pl.DataFrame,
    top_n: int = 10,
) -> alt.HConcatChart:

    data = data.with_columns(
        pl.when(pl.col("project_on") == "A")
        .then(pl.col("name_B"))
        .otherwise(pl.col("name_A"))
        .alias("projected"),
    )

    min_residue_value = int(data["residue"].min())
    max_residue_value = int(data["residue"].max())

    # -------------------------
    # Controls
    # -------------------------

    search = alt.param(
        name="search",
        value="",
        bind=alt.binding(
            input="search",
            placeholder="Search interaction...",
            name="Highlight: ",
        ),
    )

    min_residue = alt.param(
        name="min_residue",
        value=min_residue_value,
        bind=alt.binding_range(
            min=min_residue_value,
            max=max_residue_value,
            step=1,
            name="Min residue: ",
        ),
    )

    max_residue = alt.param(
        name="max_residue",
        value=max_residue_value,
        bind=alt.binding_range(
            min=min_residue_value,
            max=max_residue_value,
            step=1,
            name="Max residue: ",
        ),
    )

    min_icp = alt.param(
        name="min_icp",
        value=0.6,
        bind=alt.binding_range(
            min=0,
            max=1,
            step=0.01,
            name="Min iCP: ",
        ),
    )

    residue_filter = "datum.residue >= min_residue && datum.residue <= max_residue"

    search_match = (
        "(search === '') || (indexof(lower(datum.projected), lower(search)) >= 0)"
    )

    # -------------------------
    # Shared filtered plot data
    # -------------------------

    # First restrict to the requested residue window.
    #
    # Then compute each interaction's maximum iCP INSIDE that
    # window and keep the interaction only if that maximum
    # exceeds min_icp.
    #
    # This means we don't chop holes into the lines.
    plot_base = (
        alt.Chart(data)
        .transform_filter(residue_filter)
        .transform_joinaggregate(
            max_icp="max(icp)",
            groupby=["projected"],
        )
        .transform_filter("datum.max_icp >= min_icp")
        .encode(
            x=alt.X(
                "residue:Q",
                title="Residue #",
            ),
            y=alt.Y(
                "icp:Q",
                title="Contact probability",
                scale=alt.Scale(domain=[0, 1]),
            ),
            color=alt.Color(
                "projected:N",
                legend=None,
            ),
        )
    )

    # -------------------------
    # Lines
    # -------------------------

    lines = plot_base.mark_line(strokeWidth=1.5).encode(
        opacity=alt.condition(
            search_match,
            alt.value(1),
            alt.value(0.04),
        ),
        strokeWidth=alt.condition(
            (
                "(search !== '') && "
                "(indexof("
                "lower(datum.projected), "
                "lower(search)"
                ") >= 0)"
            ),
            alt.value(4),
            alt.value(1.5),
        ),
    )

    # -------------------------
    # High-iCP dots
    # -------------------------

    points = (
        plot_base.transform_filter("datum.icp >= min_icp")
        .mark_point(
            filled=True,
            size=45,
        )
        .encode(
            opacity=alt.condition(
                search_match,
                alt.value(1),
                alt.value(0.04),
            ),
            size=alt.condition(
                (
                    "(search !== '') && "
                    "(indexof("
                    "lower(datum.projected), "
                    "lower(search)"
                    ") >= 0)"
                ),
                alt.value(110),
                alt.value(45),
            ),
            tooltip=[
                alt.Tooltip(
                    "projected:N",
                    title="Interaction",
                ),
                alt.Tooltip(
                    "residue:Q",
                    title="Residue",
                ),
                alt.Tooltip(
                    "icp:Q",
                    title="Contact probability",
                    format=".2f",
                ),
            ],
        )
    )

    plot = (
        alt.layer(lines, points)
        .properties(
            width=850,
            height=500,
        )
        .interactive()
    )

    # -------------------------
    # Dynamic top-hits table
    # -------------------------

    # 1. Apply current residue filter.
    # 2. Within each interaction, sort residues by iCP.
    # 3. Keep its best residue.
    # 4. Apply min-iCP threshold.
    # 5. Rank interactions by that best iCP.
    hits = (
        alt.Chart(data)
        .transform_filter(residue_filter)
        .transform_window(
            best_row="row_number()",
            sort=[
                alt.SortField(
                    field="icp",
                    order="descending",
                )
            ],
            groupby=["projected"],
        )
        .transform_filter("datum.best_row == 1")
        .transform_filter("datum.icp >= min_icp")
        .transform_window(
            hit_rank="row_number()",
            sort=[
                alt.SortField(
                    field="icp",
                    order="descending",
                )
            ],
        )
        .transform_filter(f"datum.hit_rank <= {top_n}")
    )

    table_y = alt.Y(
        "hit_rank:O",
        axis=None,
        sort="ascending",
    )

    rank_column = (
        hits.mark_text(align="right")
        .encode(
            y=table_y,
            text=alt.Text(
                "hit_rank:Q",
                format=".0f",
            ),
        )
        .properties(
            width=30,
            height=500,
            title="#",
        )
    )

    name_column = (
        hits.mark_text(align="left")
        .encode(
            y=table_y,
            text="projected:N",
        )
        .properties(
            width=170,
            height=500,
            title="Interaction",
        )
    )

    residue_column = (
        hits.mark_text(align="right")
        .encode(
            y=table_y,
            text=alt.Text(
                "residue:Q",
                format=".0f",
            ),
        )
        .properties(
            width=65,
            height=500,
            title="Residue",
        )
    )

    icp_column = (
        hits.mark_text(align="right")
        .encode(
            y=table_y,
            text=alt.Text(
                "icp:Q",
                format=".3f",
            ),
        )
        .properties(
            width=60,
            height=500,
            title="iCP",
        )
    )

    dashboard = alt.hconcat(
        plot,
        rank_column,
        name_column,
        residue_column,
        icp_column,
        spacing=8,
    ).add_params(
        search,
        min_residue,
        max_residue,
        min_icp,
    )

    return dashboard


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

    hit_genes = results_table.filter(pl.col("icp").gt(0.5)).unique(["name_A", "name_B"])
    filtered_data = results_table.join(hit_genes, on=["name_A", "name_B"], how="semi")

    if len(filtered_data) == 0:
        logger.info("No hits, so nothing to plot")
    else:
        logger.info(f"Kept {len(filtered_data.unique('name_B'))} hits to plot.")
        chart = make_chart(filtered_data)
        chart.save(
            out_html,
            embed_options={
                "renderer": "canvas",
                "actions": True,
            },
        )
        logger.info("Saved plots.")

    logger.info("All done.")


if __name__ == "__main__":
    main()

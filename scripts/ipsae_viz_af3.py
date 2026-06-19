#!/usr/bin/env python3

import polars as pl
import pandas as pd
import plotly.express as px
from pathlib import Path
import argparse
import logging

logger = logging.getLogger(__name__)


def load_ipsae(ipsae_dir: Path) -> pl.DataFrame:
    files = ipsae_dir.glob("*10.txt")
    ipsae_raw = []
    for f in files:
        logger.info(f"On file {f}")
        try:
            new_df = pl.from_pandas(pd.read_csv(f, sep="\s+"))
            if len(new_df) > 0:
                ipsae_raw.append(new_df)
        except:
            logger.warning(f"Skipped {f}")
            continue

    ipsae_raw = pl.concat(ipsae_raw)

    return (
        ipsae_raw.filter(pl.col("Type") == "max")
        .with_columns(pl.col("Model").str.split("/").list.get(-2).alias("pair"))
        .sort(by="ipSAE", descending=True)
    )


def plot_metrics(ipsaes: pl.DataFrame):
    fig = px.scatter(
        ipsaes,
        title=f"{screen} screen",
        x="ipSAE",
        y="LIS",
        color="ipTM_af",
        hover_data=["pair", "ipTM_af", "pDockQ", "pDockQ2"],
        template="plotly_white",
    )
    return fig


def main():
    logger.info("Loading scores...")
    ipsae = load_ipsae(screen)
    ipsae.write_csv(out_dir / "screen_metrics.tsv", separator="\t")

    logger.info("Scores loaded. Plotting and saving...")
    fig = plot_metrics(ipsae)
    fig.write_html(out_dir / "screen_metrics.html")

    logger.info("Done.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()

    parser.add_argument("--ipsae_dir", required=True, type=Path)
    parser.add_argument("--out_dir", required=True, type=Path)

    args = parser.parse_args()
    screen = args.ipsae_dir
    out_dir = args.out_dir

    main()

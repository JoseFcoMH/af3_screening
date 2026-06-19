#!/usr/bin/env bash

#SBATCH --job-name=af3_panopticom # job name
#SBATCH --mail-type=BEGIN,END,FAIL    # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jfmh@sund.ku.dk    # email address to receive the notification
#SBATCH -c 1    # number of requested cores
#SBATCH --mem=8gb    # total requested RAM
#SBATCH --time=7-05:05:00               # max. running time of the job, format in D-HH:MM:SS
#SBATCH --output=affast_test%j.log   # Standard output and error log, '%j' gives the job ID

source /maps/projects/dan1/apps/miniconda/latest/etc/profile.d/conda.sh
conda activate random

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

CONFIG_FILE="$(realpath "$1")"
CONFIG_DIR="$(dirname "$CONFIG_FILE")"

snakemake \
  -s "$REPO_DIR/Snakefile" \
  --profile "$REPO_DIR/profiles/danhead" \
  --directory "$CONFIG_DIR" \
  --configfile "$CONFIG_FILE"

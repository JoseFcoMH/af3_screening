#!/bin/bash

#SBATCH --job-name=af3_panopticom # job name
#SBATCH --mail-type=BEGIN,END,FAIL    # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jfmh@sund.ku.dk    # email address to receive the notification
#SBATCH -c 1    # number of requested cores
#SBATCH --mem=8gb    # total requested RAM
#SBATCH --time=7-05:05:00               # max. running time of the job, format in D-HH:MM:SS
#SBATCH --output=affast_test%j.log   # Standard output and error log, '%j' gives the job ID

export PATH="$PATH:/home/vzl634/datadir/af3_screening/scripts"

module purge
module load miniconda/latest

source /maps/projects/dan1/apps/miniconda/latest/etc/profile.d/conda.sh
conda activate random

snakemake -s /home/vzl634/datadir/af3_screening/scripts/Snakefile \
  --profile default \
  --configfile $1 \
  --rerun-incomplete

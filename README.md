# AF3 screening workflow

A snakemake workflow to screen protein-protein interactions between a list of "preys" and all proteins in one or more lists or "baits". The baits may come from a .cif structure, in which case the per-residue prey contact probability will be projected on the structure.

Uses a modified [ipSAE](###ipSAE) calculation script to allow specifying an output dir and [forked Alphafast](https://github.com/JoseFcoMH/alphafast) to allow udocker (optional, may use [original](###Alphafast) if udocker is not needed)

# Dependencies
- Docker, UDocker or Singularity to run the Alphafast container
- Conda
- The workflow must be run from a conda env with:
  - Python >= 3.10
  - Snakemake
  - Biopython
  - Polars

# References
### This repo: [Publication](https://www.cell.com/molecular-cell/fulltext/S1097-2765(26)00288-1)
### Alphafast: [Repo](https://github.com/RomeroLab/alphafast)   [Paper](https://www.biorxiv.org/content/10.64898/2026.02.17.706409v1)
### ipSAE: [Repo](https://github.com/DunbrackLab/IPSAE)   [Paper](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1)

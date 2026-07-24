# AF3 screening workflow

A snakemake workflow to screen protein-protein interactions between a list of "preys" and all proteins in one or more lists or "baits". The baits may come from a .cif structure, in which case the per-residue prey contact probability will be projected on the structure.

Uses a modified [ipSAE](###ipSAE) calculation script to allow specifying an output dir and forked [Alphafast](https://github.com/JoseFcoMH/alphafast) to allow udocker (optional, may use [original](###Alphafast) if udocker is not needed)

# References

### Alphafast: [Repo](https://github.com/RomeroLab/alphafast)   [Paper](https://www.biorxiv.org/content/10.64898/2026.02.17.706409v1)
### ipSAE: [Repo](https://github.com/DunbrackLab/IPSAE)   [Paper](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1)

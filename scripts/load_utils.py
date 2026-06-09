from pathlib import Path
import requests
import logging

logger = logging.getLogger(__name__)

dimer_template = {
    "name": "dimer_prediction",
    "sequences": [
        {"protein": {"id": "A", "sequence": ""}},
        {"protein": {"id": "B", "sequence": ""}},
    ],
    "modelSeeds": [514098765],
    "dialect": "alphafold3",
    "version": 3,
}


def get_seqs_from_struc(cif_file: str) -> dict[str, str]:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    d = MMCIF2Dict(cif_file)

    ref_to_uniprot = dict(
        zip(
            d["_struct_ref.id"],
            d["_struct_ref.pdbx_db_accession"],
        )
    )

    return {
        chain.strip(): ref_to_uniprot[ref_id]
        for chains, ref_id in zip(
            d["_struct_ref_seq.pdbx_strand_id"],
            d["_struct_ref_seq.ref_id"],
        )
        for chain in chains.split(",")
    }


def fetch_uniprot_fasta(accessions):
    query = ",".join(accessions)
    logger.info(
        f"Got {len(accesions)} accesions from input list. Fetching uniprot entries..."
    )

    response = requests.get(
        "https://rest.uniprot.org/uniprotkb/stream",
        params={
            "format": "fasta",
            "query": f"accession:({query})",
        },
    )
    response.raise_for_status()

    logger.info(f"Got {len(response.text) // 2} entries back")

    return response.text


def dict_to_fasta(
    sequences: dict[str, str],
) -> str:
    seqs = ""
    for name, seq in sequences.items():
        seqs += f">{name}\n{seq}\n"
    return seqs


def fasta_to_dict(fasta_path: Path | str) -> dict[str, str]:
    fastas = Path(fasta_path).read_text().strip().splitlines()
    d = {}
    for i in range(0, len(fastas), 2):
        uniprot = fastas[i].split("|")[1]
        name = fastas[i].split("|")[2].split(" ")[0]
        seq = fastas[i + 1]
        d[f"{uniprot}_{name}"] = seq
    return d

from pathlib import Path
import re
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

protein_drug_template = {
    "name": "protein_drug_prediction",
    "sequences": [
        {"protein": {"id": "A", "sequence": ""}},
        {"ligand": {"id": "B", "smiles": ""}},
    ],
    "modelSeeds": [514098765],
    "dialect": "alphafold3",
    "version": 3,
}


UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
    r"(?:-[0-9]+)?)$"
)


def load_smi(smi_path: Path) -> dict[str, str]:
    return {
        line.strip().split(" ")[1]: line.strip().split(" ")[0]
        for line in smi_path.read_text().splitlines()
    }


def get_uniprot_ids_from_struc(cif_file: str) -> list[str]:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    d = MMCIF2Dict(cif_file)

    ref_to_uniprot = dict(
        zip(
            d["_struct_ref.id"],
            d["_struct_ref.pdbx_db_accession"],
        )
    )

    return list({ref_to_uniprot[ref_id] for ref_id in d["_struct_ref_seq.ref_id"]})


def fetch_uniprot_fasta(accessions):
    accessions = [x.strip() for x in accessions if x.strip()]

    valid = [x for x in accessions if UNIPROT_ACCESSION_RE.match(x)]
    invalid = sorted(set(accessions) - set(valid))

    if invalid:
        logger.warning(
            f"Skipping {len(invalid)} invalid UniProt accession(s): "
            f"{', '.join(invalid)}"
        )

    if not valid:
        logger.warning("No valid UniProt accessions found.")
        return ""

    batch_size = 150
    fasta = ""
    for start in range(0, len(valid), batch_size):
        valid_batch = valid[start : start + batch_size]

        query = " OR ".join(f"accession_id:{acc}" for acc in valid_batch)

        response = requests.get(
            "https://rest.uniprot.org/uniprotkb/stream",
            params={
                "format": "fasta",
                "query": query,
            },
            timeout=60,
        )
        response.raise_for_status()

        fasta += response.text

    returned = {
        line.split("|")[1]
        for line in fasta.splitlines()
        if line.startswith((">sp|", ">tr|"))
    }

    missing = sorted(set(valid) - returned)
    if missing:
        logger.warning(
            f"{len(missing)} valid UniProt accession(s) were not returned: "
            f"{', '.join(missing)}"
        )

    logger.info(f"Retrieved {len(returned)} UniProt sequence(s)")

    return fasta


def fasta_to_dict(fasta_path: Path | str) -> dict[str, str]:
    fastas = Path(fasta_path).read_text().strip().splitlines()
    d = {}
    uniprot, name = "", ""
    for line in fastas:
        if line.startswith(">"):
            uniprot = line.split("|")[1]
            name = line.split("|")[2].split(" ")[0]
            d[f"{uniprot}_{name}"] = ""
        else:
            d[f"{uniprot}_{name}"] += line.strip()
    return d

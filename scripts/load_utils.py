from contextlib import closing
from pathlib import Path
import re
import requests
import logging
from Bio import SeqIO
from functools import cache

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


@cache
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


def accession_from_id(record_id: str) -> str:
    """Convert sp|P12345|PROT_HUMAN to P12345."""
    fields = record_id.split("|")
    if len(fields) >= 3 and fields[0] in {"sp", "tr"}:
        return fields[1]
    return record_id


def index_fasta(fasta_path: Path | str, idx_path: Path | str | None = None) -> None:
    fasta = Path(fasta_path)
    if idx_path is None:
        idx_path = fasta.parent / f"{fasta.stem}.idx"
    else:
        idx_path = Path(idx_path)
    index = SeqIO.index_db(
        str(idx_path),
        str(fasta),
        "fasta",
        key_function=accession_from_id,
    )
    logger.info(
        f"Indexed {len(index)} records at {str(idx_path)} for fasta {str(fasta)}"
    )
    index.close()


def extract_fasta_records(
    index_path: str | Path,
    accessions: list[str],
) -> tuple[list[str], str]:
    """
    Return:
        missing_accessions, found_records_as_fasta_string
    """
    requested = list(dict.fromkeys(accessions))
    missing: list[str] = []
    fasta_parts: list[str] = []

    logger.info(f"Trying to retrieve {len(requested)} accessions from {index_path}")

    with closing(SeqIO.index_db(str(index_path))) as index:
        for accession in requested:
            try:
                raw_record = index.get_raw(accession)
            except KeyError:
                missing.append(accession)
                continue

            fasta_parts.append(raw_record.decode("utf-8"))

    logger.info(f"Found {len(fasta_parts)} locally, {len(missing)} missing")

    fasta_string = "".join(fasta_parts)

    return missing, fasta_string


def records_from_api(accessions: list[str]) -> str:
    batch_size = 50
    fasta = ""
    for start in range(0, len(accessions), batch_size):
        accessions_batch = accessions[start : start + batch_size]

        query = " OR ".join(f"accession_id:{acc}" for acc in accessions_batch)

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

    missing = sorted(set(accessions) - returned)
    if missing:
        logger.warning(
            f"{len(missing)} valid UniProt accession(s) were not returned: "
            f"{', '.join(missing)}"
        )

    logger.info(f"Retrieved {len(returned)} UniProt sequence(s)")

    if (len(fasta) < 2) or (fasta[0][0] != ">"):
        logger.error(f"Fasta has an inaccessions format: \n{fasta}")
        raise ValueError("Bad API response")

    return fasta


def fetch_uniprot_fasta(
    accessions: list[str], uniprot_path: None | Path | str = None
) -> str:
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

    missing_accessions = []
    fasta = ""
    if uniprot_path is not None:
        missing_accessions, local_fasta = extract_fasta_records(
            accessions=valid, index_path=uniprot_path
        )
        fasta += local_fasta
    if len(missing_accessions) > 0:
        fasta += records_from_api(accessions=missing_accessions)

    return fasta


def fasta_to_dict(fasta_path: Path | str) -> dict[str, str]:
    fastas = Path(fasta_path).read_text().strip().splitlines()
    d = {}

    if (len(fastas) < 2) or (fastas[0][0] != ">"):
        logger.error(f"Fasta at {fasta_path} has an invalid format: \n{fastas}")
        raise ValueError("Bad fasta file")

    uniprot, name = "", ""
    for line in fastas:
        if line.startswith(">"):
            uniprot = line.split("|")[1]
            name = line.split("|")[2].split(" ")[0]
            d[f"{uniprot}_{name}"] = ""
        else:
            d[f"{uniprot}_{name}"] += line.strip()
    return d

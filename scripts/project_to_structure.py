#!/usr/bin/env python3

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

from concurrent.futures import ProcessPoolExecutor
import argparse
from pathlib import Path
from collections import defaultdict
from itertools import repeat
from dataclasses import dataclass
import logging

from project_screen import DimerPrediction

logger = logging.getLogger(__name__)


def extract_sequences_from_cif(cif_file):
    cif_dict = MMCIF2Dict(str(cif_file))
    seqs = [x for x in cif_dict["_entity_poly.pdbx_seq_one_letter_code"]]
    ids = [x for x in cif_dict["_entity_poly.pdbx_strand_id"]]
    descs = [x for x in cif_dict["_entity.pdbx_description"]]
    sequences = [
        SeqRecord(Seq(seq.replace("\n", "").upper()), id=i, description=d)
        for seq, i, d in zip(seqs, ids, descs)
    ]

    return sequences


def get_best_alignment(chain: str, seq: SeqRecord):
    aligner = PairwiseAligner()
    aligner.mismatch_score = -np.inf  # Just put a gap instead
    aligner.match_score = 1
    try:
        alignment = next(aligner.align(chain, seq))
        if alignment.score > (
            len(chain) * 0.6
        ):  # Arbitrarily chosen by none other than myself
            return alignment.score, seq, alignment.aligned
        else:
            return 0.0, seq, ((), ())
    except:
        logger.warning(f"Exception occurred while trying to align to {seq.id}")
        return 0.0, seq, ((), ())


@dataclass
class AlignableDimer(DimerPrediction):
    chains: dict[str, str] | None = None
    best_score: float = 0.0
    best_seq: SeqRecord | None = None
    best_alignment = None

    def __post_init__(self):
        super().__post_init__()
        self.load_chains()
        self.project_on = {"B": "A", "A": "B"}[self.project_on]
        self.icp = self.project_icp(on=self.project_on)

    @property
    def interactor(self):
        prots = self.name.split("__")
        return prots[np.argmax(np.array(prots) != self.reference_chain)]

    @property
    def structure(self) -> Path:
        return self.prediction / f"{self.name}_model.cif"

    @property
    def match(self):
        if self.best_seq is not None:
            return self.best_seq.id, self.best_score
        return None, self.best_score

    @property
    def match_id(self):
        return self.match[0]

    def load_chains(self) -> None:
        cif_dict = MMCIF2Dict(str(self.structure))

        residues_by_entity: dict[str, list[str]] = defaultdict(list)

        for entity_id, monomer in zip(
            cif_dict["_entity_poly_seq.entity_id"],
            cif_dict["_entity_poly_seq.mon_id"],
            strict=True,
        ):
            residues_by_entity[entity_id].append(seq1(monomer, undef_code="X"))

        sequence_by_entity = {
            entity_id: "".join(residues)
            for entity_id, residues in residues_by_entity.items()
        }

        self.chains = {
            chain_id: sequence_by_entity[entity_id]
            for chain_id, entity_id in zip(
                cif_dict["_struct_asym.id"],
                cif_dict["_struct_asym.entity_id"],
                strict=True,
            )
            if entity_id in sequence_by_entity
        }

    def map_interactor(self, L: list[SeqRecord], n_cpus=1) -> None:
        assert self.chains is not None

        with ProcessPoolExecutor(max_workers=n_cpus) as executor:
            results = list(
                executor.map(
                    get_best_alignment, repeat(self.chains[self.project_on]), L
                )
            )

        for score, seq, alignment in results:
            if score > self.best_score:
                self.best_score = score
                self.best_seq = seq
                self.best_alignment = alignment

        if self.best_seq is None:
            logger.info(
                f"{self.name} chain {self.project_on} could not be aligned to anything"
            )
        else:
            logger.info(
                f"{self.name} chain {self.interactor} best aligned to {self.best_seq.id}"
            )

    def calc_attr(self) -> None:
        if self.best_seq is None or self.best_alignment is None:
            logger.info(f"No mapped attribute for {self.name}")
            return

        score_new = np.zeros(len(self.best_seq.seq))
        score_new[:] = np.nan
        for seqA_chunk, seqB_chunk in zip(
            self.best_alignment[0], self.best_alignment[1]
        ):
            score_new[seqB_chunk[0] : seqB_chunk[1]] = self.icp[
                seqA_chunk[0] : seqA_chunk[1]
            ]
        self.iCP_attr = score_new


@dataclass
class ChainMapping:
    interactor: str = ""
    best_score: float | np.float32 = 0.0
    aligned_pair: AlignableDimer | None = None


def main():

    parser = argparse.ArgumentParser(
        prog="Process AF3 inference.",
        description="Problems? Ask Jose.",
        epilog="Hope it works.",
    )
    parser.add_argument(
        "--af3_list",
        required=True,
        nargs="+",
        type=Path,
        help="Confidence jsons from predictions to align",
    )
    parser.add_argument(
        "--reference_chain", required=True, type=str, help="Chain to project"
    )
    parser.add_argument(
        "--out_file",
        required=True,
        type=Path,
        help="File where to save the output.",
    )
    parser.add_argument(
        "--cif_file", required=True, type=Path, help="Full path to .cif file to paint."
    )
    parser.add_argument(
        "--n_cpus", required=False, type=int, default=1, help="N cores to use"
    )

    in_args = parser.parse_args()

    n_cpus = in_args.n_cpus
    logger.info(f"Using {n_cpus} CPUs.")

    out_file = in_args.out_file
    cif_file = in_args.cif_file
    reference_chain = in_args.reference_chain

    refSeqs = extract_sequences_from_cif(cif_file)
    cif_name = cif_file.stem

    results = [r.parent for r in in_args.af3_list]
    results.sort()

    logger.info(f"Found {len(results)} predictions. Loading and processing...")

    af3_pairs = []
    for pair in results:
        new_pair = AlignableDimer(prediction=pair, reference_chain=reference_chain)
        new_pair.map_interactor(refSeqs, n_cpus=n_cpus)
        new_pair.calc_attr()
        af3_pairs.append(new_pair)

    assert len(af3_pairs) > 0

    logger.info(
        f"Done. Got {len(af3_pairs)} complexes, {sum(pair.best_seq is not None for pair in af3_pairs)} of which mapped to a chain in {cif_name}. Creating attribute file..."
    )

    best_dict = {
        chain: ChainMapping() for chain in set([x.match_id for x in af3_pairs])
    }
    for pair in af3_pairs:
        if pair.best_score > best_dict[pair.match_id].best_score:
            best_dict[pair.match_id] = ChainMapping(
                pair.interactor,
                pair.best_score,
                pair,
            )

    assert len(set(best_dict.keys())) == len(
        set([x.interactor for x in best_dict.values()])
    )

    best_pairs = [
        x.aligned_pair for x in best_dict.values() if x.aligned_pair is not None
    ]

    with open(out_file, "w+") as f:
        prey, bait = out_file.parents[1].name.split("__")
        f.write(f"attribute: iCP_{prey}_{bait}\nmatch mode: any\nrecipient: residues\n")
        for pair in best_pairs:
            for i, icp_val in enumerate(pair.iCP_attr):
                if not np.isnan(icp_val):
                    f.write(f"\t/{pair.match_id}:{i + 1}\t{icp_val}\n")

    logger.info(f"Wrote {str(out_file)}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

    logger.info("Mapping screen to structure...")

    main()

    logger.info("All done!")

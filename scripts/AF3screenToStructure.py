#!/usr/bin/env python3

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

import json
import os
from concurrent.futures import ProcessPoolExecutor
import argparse
from pathlib import Path
import logging

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


class af3_pair:
    def __init__(self, pred_dir=None):
        if pred_dir is not None:
            self.name = pred_dir.name
            self.pred_dir = pred_dir
            self.structure = pred_dir / f"{pred_dir.name}_model.cif"
            self.interactor = pred_dir.name.split("__")[-1]
            self.summaries = list(pred_dir.glob("*summary_confidences.json"))

        else:
            self.name = ""
            self.pred_dir = ""
            self.structure = {}
            self.interactor = ""
            self.summaries = []

        self.min_paes = []
        self.iptms = []
        self.ptms = []
        self.chains = []
        self.confidences = None
        self.iCP = None
        self.best_score = -np.inf
        self.best_seq = None
        self.best_alignment = None
        self.iCP_attr = None

    def load_chains(self):
        cif_dict = MMCIF2Dict(str(self.structure))

        three_to_one = {
            "ALA": "A",
            "CYS": "C",
            "ASP": "D",
            "GLU": "E",
            "PHE": "F",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LYS": "K",
            "LEU": "L",
            "MET": "M",
            "ASN": "N",
            "PRO": "P",
            "GLN": "Q",
            "ARG": "R",
            "SER": "S",
            "THR": "T",
            "VAL": "V",
            "TRP": "W",
            "TYR": "Y",
            "SEC": "U",
            "PYL": "O",
            "ASX": "B",
            "GLX": "Z",
            "XAA": "X",
            "C": "C",
            "U": "U",
            "T": "T",
            "G": "G",
            "A": "A",
        }

        chainID = [int(x) for x in cif_dict["_entity_poly_seq.entity_id"]]
        aas = [three_to_one[x] for x in cif_dict["_entity_poly_seq.mon_id"]]

        chainID.append("")
        aas.append("")

        i0 = chainID[0]
        currentChain = ""
        for i, aa in zip(chainID, aas):
            if i == i0:
                currentChain += aa
            else:
                self.chains.append(currentChain)
                currentChain = aa
                i0 = i

        return self

    def calc_iCP(self):
        with open(
            Path(self.pred_dir) / f"{Path(self.pred_dir).stem}_confidences.json", "r"
        ) as f:
            self.confidences = json.load(f)
            chain_len = {
                chain: self.confidences["token_chain_ids"].count(chain)
                for chain in ["A", "B"]
            }
            contact_probs = (
                np.array(self.confidences["contact_probs"])[
                    chain_len["A"] :, : chain_len["A"]
                ]
                + np.array(self.confidences["contact_probs"])[
                    : chain_len["A"], chain_len["A"] :
                ].T
            ) / 2
            self.iCP = contact_probs.max(1)

        return self

    @staticmethod
    def get_best_alignment(chain: str, seq: SeqRecord):
        aligner = PairwiseAligner()
        aligner.mismatch_score = -np.inf  # Just put a gap instead
        aligner.match_score = 1
        try:
            alignment = next(aligner.align(chain, seq))
            if alignment.score > (
                len(chain) * 0.8
            ):  # Arbitrarily chosen by none other than myself
                return alignment.score, seq, alignment.aligned
            else:
                return -np.inf, seq, ((), ())
        except:
            return -np.inf, seq, ((), ())

    @staticmethod
    def get_alignment_wrapper(args):
        chain, seq = args
        return af3_pair.get_best_alignment(chain, seq)

    def map_interactor(self, L, n_cpus=1):
        best_score = self.best_score
        best_seq = self.best_seq
        best_alignment = self.best_alignment

        args_list = [(self.chains[1], seq) for seq in L]

        with ProcessPoolExecutor(max_workers=n_cpus) as executor:
            results = list(executor.map(af3_pair.get_alignment_wrapper, args_list))

            for score, seq, alignment in results:
                if score > best_score:
                    best_score = score
                    best_seq = seq
                    best_alignment = alignment

        self.best_score = best_score
        self.best_seq = best_seq
        self.best_alignment = best_alignment

        return self

    def calc_attr(self, L=None):

        if self.best_alignment is None and L is not None:
            self.map_interactor(L)
        if self.iCP is None:
            self.calc_iCP()
        if self.best_seq is None:
            logger.info(f"No mappable attribute for {self.name}")
            return self

        score_new = np.zeros(len(self.best_seq.seq))
        score_new[:] = np.nan
        for seqA_chunk, seqB_chunk in zip(
            self.best_alignment[0], self.best_alignment[1]
        ):
            score_new[seqB_chunk[0] : seqB_chunk[1]] = self.iCP[
                seqA_chunk[0] : seqA_chunk[1]
            ]
        self.iCP_attr = score_new

        return self

    def get_attr(self):
        return self.iCP_attr

    def get_interactor(self):
        return self.interactor

    def get_match(self):
        if self.best_seq is not None:
            return self.best_seq.id, self.best_score
        return None, self.best_score

    def get_match_id(self):
        if self.best_seq is not None:
            return self.best_seq.id
        return None

    def get_match_score(self):
        return self.best_score


def main():

    n_cpus = (
        int(os.environ.get("SLURM_CPUS_ON_NODE", 0))
        or int(os.environ.get("SLURM_CPUS_PER_TASK", 0))
        or int(os.environ.get("PBS_NUM_PPN", 0))
        or os.cpu_count()  # fallback to all available CPUs and hope for the best
        or 1
    )

    logger.info(f"Using {n_cpus} CPUs.")

    parser = argparse.ArgumentParser(
        prog="Process AF3 inference.",
        description="Problems? Ask Jose.",
        epilog="Hope it works.",
    )
    parser.add_argument(
        "--af3_dir",
        required=True,
        type=Path,
        help="Directory containing the raw AF3 inference output.",
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

    in_args = parser.parse_args()

    baseDir = in_args.af3_dir
    out_file = in_args.out_file
    cif_file = in_args.cif_file

    refSeqs = extract_sequences_from_cif(cif_file)
    cif_name = cif_file.stem

    results = list(baseDir.glob("*"))
    results.sort()

    logger.info(
        f"Found {len(results)} predictions at {baseDir}. Loading and processing..."
    )
    af3_pairs = [
        af3_pair(pair).load_chains().map_interactor(refSeqs, n_cpus).calc_attr()
        for pair in results
        if (pair / f"{pair.name}_ranking_scores.csv").exists()
    ]

    assert len(af3_pairs) > 0

    logger.info(
        f"Done. Got {len(af3_pairs)} complexes, {sum(pair.best_seq is not None for pair in af3_pairs)} of which mapped to a chain in {cif_name}. Creating attribute file..."
    )

    # Messy chunks. Maybe I'll fix someday
    best_dict = {
        chain: ("", 0, None) for chain in set([x.get_match_id() for x in af3_pairs])
    }
    for pair in af3_pairs:
        match_chain_id = pair.get_match_id()
        current_match_score = pair.get_match_score()
        if current_match_score > best_dict[match_chain_id][1]:
            best_dict[match_chain_id] = (
                pair.get_interactor(),
                current_match_score,
                pair,
            )

    assert len(set(best_dict.keys())) == len(set([x[0] for x in best_dict.values()]))

    best_pairs = [x[2] for x in best_dict.values()]

    with open(out_file, "w+") as f:
        prey, bait = out_file.parents[1].name.split("__")
        f.write(f"attribute: iCP_{prey}_{bait}\nmatch mode: any\nrecipient: residues\n")
        for pair in best_pairs:
            if pair is not None:
                icp_vals = pair.get_attr()
                match_chain_id = pair.get_match_id()
                for i in range(len(icp_vals)):
                    if np.isnan(icp_vals[i]):
                        continue
                    f.write(f"\t/{match_chain_id}:{i + 1}\t{icp_vals[i]}\n")

    logger.info(f"Wrote {str(out_file)}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

    logger.info("Mapping screen to structure...")

    main()

    logger.info("All done!")

from pathlib import Path
import polars as pl

from scripts.load_utils import (
    extract_fasta_records,
    index_fasta,
    fasta_to_dict,
    fetch_uniprot_fasta,
    get_uniprot_ids_from_struc,
)

if PROTEINLIST_BAITS:

    rule make_baits_from_proteinlist:
        input:
            source=lambda wc: BAITS_DICT[wc.bait],
        output:
            accessions=f"{PAIR_DIR}/baits_accessions.txt",
        log:
            "<logs>/make_baits_from_proteinlist/<pair_name>.log",
        wildcard_constraints:
            bait="|".join(Path(b).stem for b in PROTEINLIST_BAITS),
        shell:
            r"""
            exec > >(tee {log:q}) 2>&1
            set -euo pipefail

            sort -u {input.source:q} >{output.accessions:q}
            echo "Number of unique IDs in {output.accessions:q}: $(wc -l {output.accessions:q})"
            """


if STRUCTURE_BAITS:

    rule make_baits_from_structure:
        input:
            source=lambda wc: BAITS_DICT[wc.bait],
        output:
            accessions=f"{PAIR_DIR}/baits_accessions.txt",
        log:
            "<logs>/make_baits_from_structure/<pair_name>.log",
        wildcard_constraints:
            bait="|".join(Path(b).stem for b in STRUCTURE_BAITS),
        run:
            with rule_logger(logfile=log[0], rule_name=rule) as logger:
                uni_ids = get_uniprot_ids_from_struc(input.source)
                Path(output.accessions).write_text("\n".join(uni_ids) + "\n")


rule fetch_uniprot:
    output:
        database=UNIPROT_LOCAL,
    log:
        "<logs>/fetch_uniprot.log",
    shell:
        r"""
        exec > >(tee {log:q}) 2>&1
        set -euo pipefail

        wget -qO- \
            https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
            | gunzip -c >{output.database:q}
        """


rule make_uniprot_index:
    input:
        database=ancient(UNIPROT_LOCAL),
    output:
        index=str(Path(UNIPROT_LOCAL).parent / f"{Path(UNIPROT_LOCAL).stem}.idx"),
    log:
        "<logs>/make_uniprot_index.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            index_fasta(fasta_path=input.database, idx_path=output.index)


rule make_master_fasta:
    input:
        bait_lists=expand(
            f"{PAIR_DIR}/baits_accessions.txt", prey=PREYS, bait=BAITS_NAMED
        ),
        prey_list=config["proteins_screened"],
        uniprot_db=ancient(UNIPROT_LOCAL),
        uniprot_idx=ancient(rules.make_uniprot_index.output.index),
    output:
        accessions="<work>/all_accessions.txt",
        fasta="<work>/all_chains.fasta",
        fasta_index="<work>/all_chains.idx",
    log:
        "<logs>/make_master_fasta.log",
    params:
        prey_accessions=PREYS,
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            all_accessions = set(params.prey_accessions)
            for blist in input.bait_lists:
                baccesssions = Path(blist).read_text().strip().splitlines()
                all_accessions.update(baccesssions)
            sorted_accessions = sorted(all_accessions)
            Path(output.accessions).write_text("\n".join(sorted_accessions))
            fasta = fetch_uniprot_fasta(sorted_accessions, input.uniprot_idx)
            Path(output.fasta).write_text(fasta)
            index_fasta(fasta_path=output.fasta, idx_path=output.fasta_index)


rule make_prey_fasta_per_pair:
    input:
        fasta_index=ancient(rules.make_master_fasta.output.fasta_index),
        prey_fasta=branch(
            lambda wc: exists(f"{WORK_DIR}/{wc.bait}__{wc.prey}/prey.fasta"),
            ancient(config["proteins_screened"]),
            config["proteins_screened"],
        ),
    output:
        prey_fasta="<pair_dir>/prey.fasta",
    log:
        "<logs>/make_prey_fasta_per_pair/<pair_name>.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            _, prey_fasta = extract_fasta_records(
                index_path=input.fasta_index, accessions=[wildcards.prey]
            )
            Path(output.prey_fasta).write_text(prey_fasta)


rule make_baits_fasta_per_pair:
    input:
        fasta_index=ancient(rules.make_master_fasta.output.fasta_index),
        baits_accessions="<pair_dir>/baits_accessions.txt",
    output:
        baits_fasta="<pair_dir>/baits.fasta",
    log:
        "<logs>/make_baits_fasta_per_pair/<pair_name>.log",
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            baits = Path(input.baits_accessions).read_text().strip().splitlines()
            _, baits_fasta = extract_fasta_records(
                index_path=input.fasta_index, accessions=baits
            )
            Path(output.baits_fasta).write_text(baits_fasta)


rule make_manifest:
    input:
        prey_fasta="<pair_dir>/prey.fasta",
        baits_fasta="<pair_dir>/baits.fasta",
    output:
        manifest="<pair_dir>/prediction_manifest.txt",
    log:
        "<logs>/make_manifests/<pair_name>.log",
    params:
        max_len=int(config.get("total_max_res", 2549)),
        done_preds=sorted(DB_PREDS),
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            prey_name, prey_seq = next(
                iter(fasta_to_dict(input.prey_fasta).items())
            )
            baits = fasta_to_dict(input.baits_fasta)
            predictions = []
            for bait_name, bait_seq in baits.items():
                name = "__".join(sorted([prey_name, bait_name]))
                if (
                    (total_len := len(prey_seq) + len(bait_seq)) > params.max_len
                ) and name not in params.done_preds:
                    logger.warning(
                        f"Skipped pair {name} (doesn't already exist and it's too long - {total_len} residues total, max is {params.max_len})"
                    )
                    continue
                predictions.append(name)
            logger.info(f"Wrote {len(predictions)} entries in manifest.")
            Path(output.manifest).write_text("\n".join(predictions))


checkpoint make_project_manifest:
    input:
        manifests=expand(
            f"{PAIR_DIR}/prediction_manifest.txt",
            bait=BAITS_NAMED,
            prey=PREYS,
        ),
    output:
        project_manifest="<work>/project_manifest.tsv",
    log:
        "<logs>/make_project_manifest.log",
    priority: 100
    run:
        with rule_logger(logfile=log[0], rule_name=rule) as logger:
            logger.info(f"Found {len(input.manifests)} manifests to concatenate.")
            manifest_tables = []
            for mani in input.manifests:
                bait, prey = Path(mani).parent.name.split("__")
                mani_content = Path(mani).read_text().strip().splitlines()
                mani_table = pl.DataFrame(
                    {"prediction": mani_content}
                ).with_columns(
                    pl.lit(str(mani)).alias("from_manifest"),
                    pl.lit(bait).alias("bait"),
                    pl.lit(prey).alias("prey"),
                    pl.lit(f"{bait}__{prey}").alias("pair"),
                )
                manifest_tables.append(mani_table)
            project_manifest = pl.concat(manifest_tables).sort(
                ["from_manifest", "pair", "prediction"]
            )
            project_manifest.write_csv(output.project_manifest, separator="\t")
            logger.info(
                f"Found {len(project_manifest)} requested predictions, {len(project_manifest.unique('prediction'))} of them unique."
            )

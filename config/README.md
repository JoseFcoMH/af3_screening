# Usage

Using the config file in the test directory as a base:

```yaml
reference_proteins: from_proteinlist/proteinlists.txt 
reference_structures: /path/to/structure/list.txt
proteins_screened: from_proteinlist/preys.txt

total_max_res: 4000

work_dir: from_proteinlist
scratch_dir : /tmp/af3_pipeline/test
predictions_database: /tmp/af3_pipeline/test/test_db
uniprot_path: /tmp/af3_pipeline/test/uniprot_db/uniprot_sprot.fasta

alphafast:
  script_dir: /work/repos/alphafast/scripts
  container: romerolabduke/alphafast:latest
  db_dir: /work/datadir/mmseqs2_dbs
  weights_dir: /work/datadir/
  temp_dir: /tmp/alphafast_tmp_test
  jax_compilation_cache_dir: /tmp/alphafast_jax_cache_test
  backend: udocker
  gpu_devices: "0"

modules:
  apply_defattr:
    - chimerax

  run_alphafast:
    - udocker
```

### Params

The main ones are:

  - **reference_proteins**: path to a file containing, on each line, one path to a file containing a list of uniprot accessions (also one per line).

  - **reference_structures**: path to a file containing, on each line, one path to a structure in .cif format. Sequences will be extracted from the file, and each prey will be projected to it too.

  - **proteins_screened**: path to a file containing a list of uniprot accessions. Predicitons will be the product of this list and the concatenated reference_proteins.
  
  - **total_max_res**: AF inputs whose chains surpass this threshold won't be considered for prediction.

  - **work_dir**: working directory. Outputs will be saved here.

  - **scratch_dir**: scratch directory. Temporary.

  - **predictions_database**: path to a directory containing compressed AF outputs.

  - **uniprot_path**: path where uniprot-swissprot will be saved to.

The alphafast params mirror those described in the [alphafast documentation](https://github.com/RomeroLab/alphafast).

Modules can be specified per-rule, so that the rule will load them using `module load <module>`

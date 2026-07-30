# Usage

It's probably best to specify at least the mmseqs database path and placing it in persistent storage, since it takes a long time to set up.
A persistent predictions database will avoid running the same prediction across different workflow invocations, provided they share the same database directory.

### Params

Using the config file in the test directory as a base:

```yaml
reference_proteins: proteinlists.txt 
reference_structures: structurelist.txt
proteins_screened: preys.txt

total_max_res: 4000
prediction_batch_size: 4000

out_dir: .test
scratch_dir : /tmp/af3_pipeline/test
predictions_database: /tmp/af3_pipeline/test/test_db
uniprot_path: /tmp/af3_pipeline/test/uniprot_db/uniprot_sprot.fasta

af_repo_dir: /work/repos/alphafast
af_mmseqs_db: /work/datadir/mmseqs2_dbs
af_weights_dir: /work/datadir/
```


  - **reference_proteins**: path to a file containing, on each line, one path to a file containing a list of uniprot accessions (also one per line).

  - **reference_structures**: path to a file containing, on each line, one path to a structure in .cif format. Sequences will be extracted from the file, and each prey will be projected to it too.

  - **proteins_screened**: path to a file containing a list of uniprot accessions. Predicitons will be the product of this list and the concatenated reference_proteins.
  
  - **total_max_res**: AF inputs whose chains surpass this threshold won't be considered for prediction (optional, default 4000).

  - **prediction_batch_size**: max pairs to predict at a time (optional, default 4000).

  - **out_dir**: Outputs will be saved here.

  - **scratch_dir**: scratch directory. Work will happen mostly here, so it should be local (optional, default out_dir/tmp).

  - **predictions_database**: path to a directory containing compressed AF outputs (optional).

  - **uniprot_path**: path where uniprot-swissprot will be saved to (optional).

  - **af_repo_dir**: alphafast repo dir (optional).

  - **af_mmseqs_db**: path to mmseqs db (optional).

  - **af_weights_dir**: path to directory containing Alphafold 3 weights.

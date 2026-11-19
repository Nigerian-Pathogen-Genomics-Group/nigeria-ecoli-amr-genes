# Running Bactopia on Seqera Cloud (issue #8)

Bactopia was run from the Seqera Cloud Launchpad, so there is no single command line. The settings that were used
are in these files:

| File | Purpose |
|---|---|
| `download_assemblies.sh` | downloads the 125 assemblies from NCBI (list: `data/bactopia_accessions_125.txt`) |
| `make_sample_table.sh` | builds the Bactopia sample table (`runtype = assembly`) |
| `pre_run_script.sh` | copies the table to the head node and fills the dataset cache |
| `seqera_nextflow.config` | extra Nextflow config (plugin, docker registry, params) |
| `../../results/bactopia/run_parameters.json` | pipeline parameters |

Launchpad settings: pipeline `https://github.com/bactopia/bactopia`, revision `v4.1.0`, Nextflow 26.10,
Nextflow syntax parser v2, config profile `docker`, Seqera Compute (AWS Batch, eu-west-2, Fusion, Wave).

Illustrative command line for readers who do not use Seqera (not the one used, and not tested):

```bash
nextflow run bactopia/bactopia -r v4.1.0 -profile docker \
    --samples samples_all120.txt --outdir bactopia-all120 \
    --ask_merlin --ectyper_verify false \
    --amrfinderplus_organism Escherichia --mlst_scheme escherichia
```

Lessons learned (for anyone repeating this on Seqera Cloud):
- The sample table must be a local file on the head node (an `s3://` path failed).
- Set `datasets_cache` to a bucket path and fill it first (see `pre_run_script.sh`); downloading in-pipeline timed out.
- The `docker` profile is needed; the default profile uses conda.
- `--ectyper_verify` stops the whole run on a non-*E. coli* genome (one *Shigella* did). Use `false` for datasets from NCBI groups that mix *E. coli* and *Shigella*.
- Resuming a failed run (same session) reuses finished work, so the failed run's cost is not repeated.

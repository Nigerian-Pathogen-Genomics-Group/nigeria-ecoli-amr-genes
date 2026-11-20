# Reproducibility Record: Bactopia Analysis of 125 Nigerian E. coli and Shigella Isolates

## 1. Overview

This document describes the computational environment, software versions, databases, input genomes, parameters and execution procedures used for the Bactopia analysis of 125 Nigerian Escherichia coli and Shigella isolates.

The analysis was performed using Bactopia v4.1.0 through Seqera Cloud, with Nextflow DSL2 and AWS Batch infrastructure.

The analysis consisted of two runs:

- Initial test run: 5 isolates. All tasks completed successfully.
- Main analysis: 120 isolates. The initial execution encountered an ECTyper verification error for a Shigella genome. The workflow was resumed with ECTyper verification disabled, and all tasks subsequently completed successfully.

Together, the two runs covered 125 isolates.

## 2. Software versions

The following versions were recorded in the analysis documentation:

| Software | Version |
|---|---|
| Bactopia | 4.1.0 |
| Nextflow | 26.10 |
| nf-bactopia plugin | 2.1.7 |
| AMRFinderPlus | 4.2.7 |
| mlst | 2.35.0 |
| ECTyper | 2.0.0 |
| PathoType | v2024-10-03 |
| ClermontTyping | 24.02 |
| STECFinder | 1.1.2 |
| ShigaTyper | 2.0.5 |
| ShigaPass | 1.5.0 |

Note: These are the versions reported in the existing analysis record. Exact patch versions and container image digests remain to be verified against the original execution logs.

## 3. Computational environment

The analysis was executed using:

| Component | Recorded configuration |
|---|---|
| Workflow platform | Seqera Cloud |
| Compute platform | Seqera Compute |
| Cloud provider | AWS |
| AWS region | eu-west-2 |
| Executor | AWS Batch |
| Workflow language | Nextflow DSL2 |
| Container engine | Docker |
| Container registry | quay.io |
| Nextflow configuration profile | docker |
| AWS Batch queue | general |
| Bactopia plugin | nf-bactopia 2.1.7 |
| AWS plugin | nf-amazon |

Additional Nextflow configuration:

```groovy
plugins {
    id 'nf-amazon'
    id 'nf-bactopia@2.1.7'
}

docker {
    enabled = true
    registry = 'quay.io'
}

params {
    queue = 'general'
    disable_scratch = false
    max_memory = '128.GB'
    max_time = '240.h'
}
```

The maximum memory and execution time above are configured limits, not measured resource consumption.

The following information remains to be recovered from Seqera:

- AWS instance types.
- Actual CPU allocation.
- Peak memory usage.
- Operating system.
- Docker version.
- Actual execution resource consumption.

## 4. Reference databases

The pre-run script used to prepare the Bactopia database cache identifies the following database resources:

| Database/resource | Recorded identifier |
|---|---|
| Bactopia datasets | v4.1.0 |
| AMRFinderPlus database | Bactopia v4.1.0 dataset |
| MLST database | Bactopia v4.1.0 dataset |
| GTDB | rs207 |
| Mash reference database | RefSeq 88 |
| Mash k-mer database | k21 |

The following database files were explicitly referenced:

```text
amrfinderplus.tar.gz
mlst.tar.gz
gtdb-rs207.genomic-reps.dna.k31.lca.json.gz
mash-refseq88.k21.msh.xz
```

The original pre-run script retrieved database resources from:

```text
https://datasets.bactopia.com/datasets
```

Database files were cached in an S3 bucket to avoid repeated downloads during workflow execution.

The original analyst reported that downloading datasets during pipeline execution was too slow on Fusion and exceeded the 8-hour execution limit.

Exact database download dates, checksums and any additional database resources remain to be verified.

## 5. Input genomes and provenance

The analysis included 125 NCBI Pathogen Detection assemblies representing Nigerian E. coli and Shigella isolates.

The input assemblies were provided as FASTA files and supplied to Bactopia through sample tables.

Two sample tables were used:

- `samples_test5.txt`: 5-isolate test run.
- `samples_all120.txt`: 120-isolate main run.

The recorded main-run parameter was:

```json
{
  "samples": "samples_all120.txt"
}
```

The complete accession list and assembly versions are maintained separately in the project data files.

The following provenance information remains to be completed:

- Full accession list with assembly versions.
- Assembly download dates.
- FASTA checksums.
- Procedure for handling updated or withdrawn NCBI assemblies.

## 6. Pipeline parameters and execution configuration

The recorded main-run parameters were:

```json
{
  "samples": "samples_all120.txt",
  "outdir": "s3://<your-bucket>/bactopia-all120",
  "ask_merlin": true,
  "ectyper_verify": false,
  "amrfinderplus_organism": "Escherichia",
  "mlst_scheme": "escherichia",
  "mykrobe_species": "tb",
  "datasets_cache": "s3://<your-bucket>/bactopia-datasets"
}
```

The pipeline repository was:

https://github.com/bactopia/bactopia

Pipeline revision:

```text
v4.1.0
```

Configuration profile:

```text
docker
```

Nextflow syntax parser:

```text
v2
```

The initial 5-isolate test used the same principal parameters, except:

```text
ectyper_verify=true
samples=samples_test5.txt
```

### Pre-run preparation

Before the main analysis, the samplesheet was copied from cloud storage into the local execution environment:

```bash
aws s3 cp s3://<your-bucket>/inputs/all120/samples_all120.txt \
    /tmp/samples_all120.txt
```

The database cache was prepared using a separate pre-run script.

The script checked whether the required database files already existed in S3 and downloaded missing resources.

The actual S3 bucket names have been replaced with placeholders in this document.

### ECTyper verification

During the initial execution of the 120-isolate run, ECTyper's `--verify` option rejected a Shigella genome:

```text
GCA_013391325
```

The workflow was subsequently resumed with:

```text
ectyper_verify=false
```

The remaining workflow tasks completed successfully.

Consequently, the ECTyper verification check was disabled for the main 120-isolate run.

This difference between the test and main runs must be retained when reproducing the analysis.

## 7. Randomness and time-dependent steps

The following factors may affect reproducibility:

1. NCBI assemblies may be updated or withdrawn. The original assembly version must therefore be recorded and used when reproducing the analysis.
2. Bactopia databases may change between releases. The original database versions should be retained.
3. Database downloads were cached in S3 to avoid repeated downloads and changes in execution behaviour.
4. The initial test and main analysis used different ECTyper verification settings.
5. The original S3 output and cache locations must be recorded for the complete execution history.

No random seeds or explicitly stochastic analysis steps have yet been identified in the available execution records.

## 8. Outputs

The original analysis produced the following outputs:

| File/directory | Description |
|---|---|
| `raw/test5/` | Original merged Bactopia results from the 5-isolate test |
| `raw/all120/` | Original merged Bactopia results from the 120-isolate run |
| `bactopia_summary_125.tsv` | Combined isolate-level summary |
| `amr_gene_matrix_125.tsv` | AMR gene and mutation presence matrix |
| `bactopia_vs_ncbi_125.tsv` | Comparison of Bactopia AMR calls with NCBI AMR genotypes |
| `run_parameters.json` | Recorded main-run parameters |

The original merged Bactopia result tables were retained without modification.

The combined summary and comparison tables were generated using project-specific Python scripts.

## 9. Outstanding information

The following information must be obtained from the original Seqera Cloud execution records before this document can be considered complete:

- Exact Nextflow version, including patch version.
- Exact software versions verified from execution logs.
- Docker version and container image digests.
- AWS instance types and CPU allocation.
- Actual memory usage and execution times.
- Operating system details.
- Database download dates and checksums.
- Complete input accession list with assembly versions.
- Input FASTA checksums.
- Exact original launch command and resolved configuration.
- Random seeds or other time-dependent execution settings, if applicable.

These details will be added after verification against the original execution records.

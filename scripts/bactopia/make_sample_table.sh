#!/usr/bin/env bash
# Build the Bactopia sample table for assemblies that were uploaded to cloud storage.
# Usage: scripts/bactopia/make_sample_table.sh inputs/fasta s3://<your-bucket>/inputs/all120 > samples.txt
# The columns are those of `bactopia prepare` for assemblies; genome_size is a nominal 5.15 Mb for E. coli.
set -euo pipefail
FASTA_DIR="$1"; REMOTE="$2"
printf 'sample\truntype\tgenome_size\tspecies\tr1\tr2\tse\tont\tassembly\n'
for f in "$FASTA_DIR"/*.fna.gz; do
  a=$(basename "$f" .fna.gz)
  printf '%s\tassembly\t5150000\tEscherichia coli\t\t\t\t\t%s/%s.fna.gz\n' "$a" "$REMOTE" "$a"
done

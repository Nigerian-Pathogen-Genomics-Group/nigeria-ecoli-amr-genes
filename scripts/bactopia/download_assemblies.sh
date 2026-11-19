#!/usr/bin/env bash
# Download NCBI assemblies (FASTA) for a list of accessions (one per line, e.g. GCA_060580055.1).
# Needs the NCBI "datasets" command-line tool (conda install -c conda-forge ncbi-datasets-cli).
# Usage: scripts/bactopia/download_assemblies.sh data/bactopia_accessions_125.txt inputs/fasta
set -euo pipefail
LIST="$1"; OUT="${2:-fasta}"
mkdir -p "$OUT"
while read -r acc; do
  [[ -z "$acc" ]] && continue
  base="${acc%%.*}"                      # drop the version for the file name
  [[ -s "$OUT/$base.fna.gz" ]] && continue
  tmp=$(mktemp -d)
  datasets download genome accession "$acc" --include genome --filename "$tmp/x.zip"
  unzip -q "$tmp/x.zip" -d "$tmp"
  gzip -c "$tmp"/ncbi_dataset/data/"$acc"/*.fna > "$OUT/$base.fna.gz"
  echo "$acc" >> "$OUT/accession_versions.txt"   # record the exact version downloaded
  rm -rf "$tmp"
done < "$LIST"

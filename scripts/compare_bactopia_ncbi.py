#!/usr/bin/env python3
"""Compare Bactopia (AMRFinderPlus) resistance calls with the NCBI Pathogen Detection calls.

Usage:
    python3 scripts/compare_bactopia_ncbi.py \
        --runs results/bactopia/raw/test5 results/bactopia/raw/all120 \
        --ncbi data/Ecoli_AMR_analysis_125.tsv \
        --out results/bactopia/bactopia_vs_ncbi_125.tsv

What is compared: the set of element symbols per isolate.
  * NCBI: the "AMR genotypes" column, entries of the form gene=STATUS; STATUS is dropped.
  * Bactopia: amrfinderplus.tsv rows with Type == AMR, column "Element symbol".
Repeat hits of the same gene in one isolate count once. Point mutations (for example
gyrA_S83L) are compared by name like any other element. Note: only gene/mutation names
are compared, not coverage or identity values. NCBI also uses AMRFinderPlus, so a match
means our run reproduces NCBI's calls; it is not an independent validation.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--ncbi", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    bact = defaultdict(set)
    for run in a.runs:
        for r in read_tsv(Path(run) / "amrfinderplus.tsv"):
            if r["Type"] == "AMR":
                bact[r["Name"]].add(r["Element symbol"])

    ncbi = {}
    for r in read_tsv(a.ncbi):
        acc = r["Assembly"].split(".")[0]
        ncbi[acc] = {g.split("=")[0] for g in r["AMR genotypes"].split(",") if g}

    isolates = sorted(set(ncbi) | set(bact))
    n_same = 0
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "n_NCBI", "n_Bactopia", "n_shared", "identical",
                    "only_in_NCBI", "only_in_Bactopia"])
        for acc in isolates:
            n, b = ncbi.get(acc, set()), bact.get(acc, set())
            same = n == b
            n_same += same
            w.writerow([acc, len(n), len(b), len(n & b), "yes" if same else "no",
                        ";".join(sorted(n - b)), ";".join(sorted(b - n))])
    total_n = sum(len(v) for v in ncbi.values())
    total_b = sum(len(bact.get(k, ())) for k in ncbi)
    print(f"Isolates compared: {len(isolates)}")
    print(f"Identical gene sets: {n_same} of {len(isolates)}")
    print(f"Total calls NCBI: {total_n}; Bactopia: {total_b}")
    print(f"Wrote {a.out}")


if __name__ == "__main__":
    main()

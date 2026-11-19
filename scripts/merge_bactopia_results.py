#!/usr/bin/env python3
"""Merge Bactopia merged-results tables from one or more runs into one table per isolate.

Usage:
    python3 scripts/merge_bactopia_results.py \
        --runs results/bactopia/raw/test5 results/bactopia/raw/all120 \
        --ncbi data/Ecoli_AMR_analysis_125.tsv \
        --outdir results/bactopia

Inputs (per run directory, original Bactopia file names, unmodified):
    amrfinderplus.tsv mlst.tsv ectyper.tsv clermontyping.tsv stecfinder.tsv
    shigapass.tsv assembly-scan.tsv  (shigatyper.tsv is optional)

Outputs (in --outdir):
    bactopia_summary_125.tsv   one row per isolate
    amr_gene_matrix_125.tsv    isolates x AMR genes, 1 = present (Type == AMR only)

Rules:
  * Isolates are matched on the assembly accession without the version suffix
    (GCA_002284675.1 -> GCA_002284675).
  * AMR genes: only rows with Type == AMR; each element symbol counts once per
    isolate (repeat hits of the same gene are collapsed). STRESS and VIRULENCE
    rows are ignored.
  * The script stops if the tools do not cover exactly the same isolates.
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(runs, name, key, strip=""):
    """Read `name` from every run directory and index rows by isolate."""
    out = {}
    found = False
    for run in runs:
        p = Path(run) / name
        if not p.exists():
            continue
        found = True
        for r in read_tsv(p):
            k = r[key]
            if strip and k.endswith(strip):
                k = k[: -len(strip)]
            if k in out:
                sys.exit(f"ERROR: {k} appears twice in {name}")
            out[k] = r
    return out if found else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--ncbi", required=True, help="NCBI Pathogen Detection table (125 isolates)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    mlst = load(a.runs, "mlst.tsv", "FILE")
    clermont = load(a.runs, "clermontyping.tsv", "sample", strip=".fna")
    stec = load(a.runs, "stecfinder.tsv", "Sample")
    ect = load(a.runs, "ectyper.tsv", "Name")
    asm = load(a.runs, "assembly-scan.tsv", "sample")
    spass = load(a.runs, "shigapass.tsv", "Name")
    styper = load(a.runs, "shigatyper.tsv", "sample") or {}

    amr = defaultdict(set)
    amr_names = set()
    for run in a.runs:
        for r in read_tsv(Path(run) / "amrfinderplus.tsv"):
            amr_names.add(r["Name"])
            if r["Type"] == "AMR":
                amr[r["Name"]].add(r["Element symbol"])

    ncbi = {r["Assembly"].split(".")[0]: r for r in read_tsv(a.ncbi)}

    sets = {"mlst": set(mlst), "clermontyping": set(clermont), "stecfinder": set(stec),
            "ectyper": set(ect), "assembly-scan": set(asm), "shigapass": set(spass),
            "amrfinderplus": amr_names, "ncbi": set(ncbi)}
    ref = sets["ncbi"]
    bad = {k: (sorted(ref - v), sorted(v - ref)) for k, v in sets.items() if v != ref}
    if bad:
        sys.exit(f"ERROR: tools do not cover the same isolates (missing, extra): {bad}")
    print(f"OK: all tools cover the same {len(ref)} isolates")

    cols = ["accession", "assembly_version", "biosample", "location", "isolation_source_ncbi",
            "taxon", "ST", "mlst_status", "ST_shigapass", "phylogroup", "phylogroup_mash",
            "O_type", "H_type", "serotype_ectyper", "stecfinder_cluster", "stx", "ipaH",
            "shigapass", "shigatyper", "ipaB", "n_AMR_genes",
            "total_length", "n_contigs", "n50"]
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for acc in sorted(ref):
        n, e, c, s, m, q = ncbi[acc], ect[acc], clermont[acc], stec[acc], mlst[acc], asm[acc]
        sp, ty = spass[acc], styper.get(acc)
        rows.append({
            "accession": acc, "assembly_version": n["Assembly"], "biosample": n["BioSample"],
            "location": n["Location"], "isolation_source_ncbi": n["Isolation source"],
            "taxon": e["Species"], "ST": m["ST"], "mlst_status": m["STATUS"],
            "ST_shigapass": sp["MLST"], "phylogroup": c["phylogroup"], "phylogroup_mash": c["mash_group"],
            "O_type": e["O-type"], "H_type": e["H-type"], "serotype_ectyper": e["Serotype"],
            "stecfinder_cluster": s["Cluster"], "stx": s["stx type"], "ipaH": s["ipaH presence"],
            "shigapass": sp["Predicted_Serotype"],
            "shigatyper": ty["prediction"] if ty else "NA", "ipaB": ty["ipaB"] if ty else "NA",
            "n_AMR_genes": len(amr[acc]), "total_length": q["total_contig_length"],
            "n_contigs": q["total_contig"], "n50": q["n50_contig_length"],
        })
    with open(outdir / "bactopia_summary_125.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    genes = sorted({g for v in amr.values() for g in v})
    with open(outdir / "amr_gene_matrix_125.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession"] + genes)
        for acc in sorted(ref):
            w.writerow([acc] + [int(g in amr[acc]) for g in genes])
    print(f"Wrote {len(rows)} isolates, {len(genes)} AMR genes/mutations to {outdir}")


if __name__ == "__main__":
    main()

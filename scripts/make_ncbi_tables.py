#!/usr/bin/env python3
"""Build the NCBI-genotype tables for the poster (issues #6, #11, #12, #14).

Run from the repository root:  python3 scripts/make_ncbi_tables.py
Reads:   data/Ecoli_AMR_analysis_125.tsv, data/source_groups_125.tsv
Writes:  results/ncbi_tables/*.tsv

Rules:
- AMR genotypes are the ones NCBI assigned (AMRFinderPlus); nothing is re-run here.
- A gene counts as present if its NCBI status is COMPLETE; a point mutation (POINT) counts as present.
  PARTIAL, PARTIAL_END_OF_CONTIG and MISTRANSLATION are NOT counted (they stay visible in the status matrix).
- Source groups come from data/source_groups_125.tsv. No statistical comparison is made between groups.
"""
import csv
import collections
import os

ANALYSIS = "data/Ecoli_AMR_analysis_125.tsv"
GROUPS = "data/source_groups_125.tsv"
OUT = "results/ncbi_tables"
ORDER = ["Human", "Animal/food-chain", "Environment"]
PRESENT = ("COMPLETE", "POINT")

os.makedirs(OUT, exist_ok=True)
rows = list(csv.DictReader(open(ANALYSIS), delimiter="\t"))
grp = {r["Assembly"]: r["Source group"] for r in csv.DictReader(open(GROUPS), delimiter="\t")}

iso = []
for r in rows:
    genes = {}
    for token in r["AMR genotypes"].split(","):
        token = token.strip()
        if token:
            name, status = token.rsplit("=", 1)
            genes[name] = status
    iso.append(dict(acc=r["Assembly"], bs=r["BioSample"], src=r["Isolation source"],
                    loc=r["Location"], grp=grp[r["Assembly"]], genes=genes))
assert len(iso) == 125 and all(i["grp"] in ORDER for i in iso)
iso.sort(key=lambda i: (ORDER.index(i["grp"]), i["acc"]))


def present(i, gene):
    return i["genes"].get(gene) in PRESENT


def write(name, header, lines):
    with open(os.path.join(OUT, name), "w") as f:
        f.write("\t".join(header) + "\n")
        for line in lines:
            f.write("\t".join(str(x) for x in line) + "\n")


# Issue #6: isolates by source group, isolation source and location
def clean_source(s):
    return "Human stool samples" if s.startswith("Human stool") else s

by_source = collections.Counter((i["grp"], clean_source(i["src"])) for i in iso)
write("table_isolates_by_source_125.tsv", ["Source group", "Isolation source", "Isolates"],
      [(g, s, n) for g in ORDER for (gg, s), n in sorted(by_source.items(), key=lambda x: -x[1]) if gg == g])
by_loc = collections.Counter((i["grp"], i["loc"]) for i in iso)
write("table_source_by_location_125.tsv", ["Source group", "Location", "Isolates"],
      [(g, l, n) for g in ORDER for (gg, l), n in sorted(by_loc.items()) if gg == g])

# Issue #11: isolates by genes (status matrix and 0/1 matrix)
genes = sorted({g for i in iso for g in i["genes"]})
head = ["Assembly", "BioSample", "Source group", "Isolation source"] + genes
write("isolates_by_genes_status_125.tsv", head,
      [[i["acc"], i["bs"], i["grp"], i["src"]] + [i["genes"].get(g, "") for g in genes] for i in iso])
write("isolates_by_genes_binary_125.tsv", head,
      [[i["acc"], i["bs"], i["grp"], i["src"]] + [int(present(i, g)) for g in genes] for i in iso])

# Issue #12: frequency by source group
n = {o: sum(1 for i in iso if i["grp"] == o) for o in ORDER}
freq = []
for g in genes:
    cnt = {o: sum(1 for i in iso if i["grp"] == o and present(i, g)) for o in ORDER}
    freq.append((g, sum(cnt.values()), cnt))
freq.sort(key=lambda x: (-x[1], x[0]))
header = ["Gene or mutation", "Type", "Total n", "Total %"]
for o in ORDER:
    header += [f"{o} n (of {n[o]})", f"{o} %"]
lines = []
for g, tot, cnt in freq:
    typ = "point mutation" if any(i["genes"].get(g) == "POINT" for i in iso) else "gene"
    line = [g, typ, tot, f"{100 * tot / 125:.1f}"]
    for o in ORDER:
        line += [cnt[o], f"{100 * cnt[o] / n[o]:.1f}"]
    lines.append(line)
write("gene_frequency_by_source_125.tsv", header, lines)

# Issue #14: isolates with mcr or blaCTX-M genes
lines = []
for i in iso:
    mcr = [f"{g}({s})" for g, s in i["genes"].items() if g.startswith("mcr")]
    ctx = [f"{g}({s})" for g, s in i["genes"].items() if g.startswith("blaCTX-M")]
    if mcr or ctx:
        lines.append([i["acc"], i["bs"], i["grp"], i["src"], i["loc"], ";".join(mcr), ";".join(ctx)])
write("mcr_and_blaCTX-M_isolates_125.tsv",
      ["Assembly", "BioSample", "Source group", "Isolation source", "Location", "mcr genes", "blaCTX-M genes"], lines)
print("done:", len(genes), "genes/mutations;", len(lines), "isolates with mcr or blaCTX-M")

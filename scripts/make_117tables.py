#!/usr/bin/env python3
"""Issue #13: the #11-style tables for the 117 E. coli (new files; the 125-isolate tables are not changed).

Run from the repo root:
    PYTHONPATH=scripts python3 scripts/make_tables_117.py
Writes into results/tables_117/:
    isolates_by_genes_binary_117.tsv   1 = element detected by AMRFinderPlus (Bactopia run), 0 = not detected
    isolates_by_genes_status_117.tsv   the NCBI call for each element: COMPLETE, POINT, PARTIAL, PARTIAL_END_OF_CONTIG, MISTRANSLATION
    isolates_by_source_117.tsv         isolates per source group and isolation source
    isolates_by_source_and_bioproject_117.tsv
    source_by_location_117.tsv
    mcr_and_blaCTX-M_isolates_117.tsv
The gene-frequency table (#12) is written by scripts/fig3_frequency.py (results/gene_frequency_by_source_117.tsv).

Rule: an element counts as present if AMRFinderPlus reported it (Type == AMR), whatever its status. This is the same
rule as the figures. NOTE: the 125-isolate tables in results/ncbi_tables count only COMPLETE and POINT calls, so numbers
for a few elements differ between those tables and these. The status table shows which calls are PARTIAL etc.
"""
import csv
from pathlib import Path
import pandas as pd
from fig_common import load, GROUP_ORDER

OUT = Path("results/tables_117")
OUT.mkdir(parents=True, exist_ok=True)
meta, X = load()

ncbi = pd.read_csv("data/Ecoli_AMR_analysis_125.tsv", sep="\t")
ncbi["acc"] = ncbi["Assembly"].str.split(".").str[0]
ncbi = ncbi.set_index("acc")
meta["assembly"] = ncbi.loc[meta.index, "Assembly"]
meta["isolation_source"] = ncbi.loc[meta.index, "Isolation source"].map(
    lambda s: "Human stool samples" if s.startswith("Human stool") else s)
meta["location"] = ncbi.loc[meta.index, "Location"]
meta["biosample_id"] = meta["biosample"]
meta["group"] = pd.Categorical(meta["group"], GROUP_ORDER, ordered=True)
meta = meta.sort_values(["group", "assembly"])
X = X.loc[meta.index]
genes = list(X.columns)

status = {}
for acc in meta.index:
    d = {}
    for tok in ncbi.loc[acc, "AMR genotypes"].split(","):
        if tok.strip():
            name, st = tok.strip().rsplit("=", 1)
            d[name] = st
    status[acc] = d

base = pd.DataFrame({"Assembly": meta["assembly"], "BioSample": meta["biosample_id"], "Source group": meta["group"].astype(str),
                     "BioProject": meta["bioproject"], "Isolation source": meta["isolation_source"]})
b = base.join(X.astype(int))
b.to_csv(OUT / "isolates_by_genes_binary_117.tsv", sep="\t", index=False)
st = base.copy()
for g in genes:
    st[g] = [status[a].get(g, "") for a in meta.index]
st.to_csv(OUT / "isolates_by_genes_status_117.tsv", sep="\t", index=False)

(meta.groupby(["group", "isolation_source"], observed=True).size().rename("Isolates").reset_index()
     .sort_values(["group", "Isolates"], ascending=[True, False])
     .rename(columns={"group": "Source group", "isolation_source": "Isolation source"})
     .to_csv(OUT / "isolates_by_source_117.tsv", sep="\t", index=False))
(meta.groupby(["group", "bioproject"], observed=True).size().rename("Isolates").reset_index()
     .sort_values(["group", "Isolates"], ascending=[True, False])
     .rename(columns={"group": "Source group", "bioproject": "BioProject"})
     .to_csv(OUT / "isolates_by_source_and_bioproject_117.tsv", sep="\t", index=False))
(meta.groupby(["group", "location"], observed=True).size().rename("Isolates").reset_index()
     .rename(columns={"group": "Source group", "location": "Location"})
     .to_csv(OUT / "source_by_location_117.tsv", sep="\t", index=False))

rows = []
for acc, r in meta.iterrows():
    mcr = [f"{g}({status[acc].get(g, 'detected')})" for g in genes if g.startswith("mcr") and X.loc[acc, g]]
    ctx = [f"{g}({status[acc].get(g, 'detected')})" for g in genes if g.startswith("blaCTX-M") and X.loc[acc, g]]
    if mcr or ctx:
        rows.append([r["assembly"], r["biosample_id"], str(r["group"]), r["bioproject"], r["isolation_source"], r["location"],
                     ";".join(mcr), ";".join(ctx)])
pd.DataFrame(rows, columns=["Assembly", "BioSample", "Source group", "BioProject", "Isolation source", "Location",
                            "mcr genes", "blaCTX-M genes"]).to_csv(OUT / "mcr_and_blaCTX-M_isolates_117.tsv", sep="\t", index=False)

n_calls = int(X.values.sum())
n_partial = sum(1 for a in meta.index for g, s in status[a].items() if g in genes and s not in ("COMPLETE", "POINT"))
print(f"{len(meta)} isolates, {len(genes)} elements, {n_calls} calls; {n_partial} calls are not COMPLETE/POINT (partial etc.)")
print("isolates with blaCTX-M:", len(rows), "; with mcr:", sum(1 for r in rows if r[6]))
print("wrote", len(list(OUT.glob('*.tsv'))), "files in", OUT)
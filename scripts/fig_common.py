#!/usr/bin/env python3
"""Shared data loading for the issue #13 figures (117 confirmed E. coli).

Run the figure scripts from the repo root, e.g.  python3 scripts/fig1_heatmap.py
Package versions used when this was written: python 3.10+, numpy, pandas, scipy, matplotlib
(print them with:  python3 -c "import numpy,pandas,scipy,matplotlib as m;print(numpy.__version__,pandas.__version__,scipy.__version__,m.__version__)")
"""
import pandas as pd

SUMMARY = "results/bactopia/bactopia_summary_125.tsv"
MATRIX = "results/bactopia/amr_gene_matrix_125.tsv"
GROUPS = "data/source_groups_125.tsv"

# BioSample prefix (first 8 characters) -> BioProject. Checked on the NCBI BioSample pages (issue #13).
# SAMN6270 samples are linked to PRJNA1071873 and PRJNA514245; the legend shows PRJNA1071873.
BIOPROJECT = {
    "SAMN0719": "PRJNA389301",
    "SAMN1105": "PRJNA231221",
    "SAMN1185": "PRJNA544810",
    "SAMN1648": "PRJNA670054",
    "SAMN5288": "PRJNA1348097",
    "SAMN6076": "PRJNA1477568",
    "SAMN6184": "PRJNA1498538",
    "SAMN6270": "PRJNA1071873",
}

# Okabe-Ito colours, checked with the palette validator (colour-blind safe).
GROUP_ORDER = ["Human", "Animal/food-chain", "Environment"]
GROUP_COLOURS = {"Human": "#0072B2", "Animal/food-chain": "#E69F00", "Environment": "#009E73"}


def load():
    """Return (meta, X): meta has one row per isolate (index = accession GCA_xxxxxxxxx);
    X is the 0/1 matrix (isolates x AMR elements) for the 117 E. coli, all-zero columns dropped."""
    summ = pd.read_csv(SUMMARY, sep="\t")
    grp = pd.read_csv(GROUPS, sep="\t")
    grp["accession"] = grp["Assembly"].str.split(".").str[0]
    grp = grp.rename(columns={"Source group": "group"})
    meta = summ.merge(grp[["accession", "group"]], on="accession", how="left")
    assert meta["group"].notna().all(), "isolate without a source group"
    meta = meta[meta["taxon"] == "Escherichia coli"].copy()
    meta["bioproject"] = meta["biosample"].str[:8].map(BIOPROJECT)
    assert meta["bioproject"].notna().all(), "BioSample prefix not in BIOPROJECT"
    meta = meta.set_index("accession")
    X = pd.read_csv(MATRIX, sep="\t").set_index("accession").loc[meta.index]
    X = X.loc[:, X.sum() > 0]
    assert len(meta) == 117, f"expected 117 isolates, got {len(meta)}"
    return meta, X


def bioproject_order(meta):
    """BioProjects ordered by number of isolates (largest first)."""
    return list(meta["bioproject"].value_counts().index)


def bioproject_colours(order):
    """Seven categorical colours (reference palette, order validated: adjacent-pair CVD and normal-vision checks pass).
    Panels using them also use a different marker shape per BioProject and a legend, so colour is never the only cue."""
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    return {bp: pal[i] for i, bp in enumerate(order)}


if __name__ == "__main__":
    meta, X = load()
    print("isolates:", len(meta), " elements:", X.shape[1])
    print(meta["group"].value_counts().to_string())
    print(pd.crosstab(meta["bioproject"], meta["group"]).to_string())
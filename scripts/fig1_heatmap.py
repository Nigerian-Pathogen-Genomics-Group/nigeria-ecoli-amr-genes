#!/usr/bin/env python3
"""Issue #13, figure 1: heatmap of AMR elements in the 117 E. coli (descriptive only).

Run from the repo root:  python3 scripts/fig1_heatmap.py
Writes figures/fig1_heatmap.pdf and figures/fig1_heatmap.png
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from fig_common import load, bioproject_order, bioproject_colours, GROUP_ORDER, GROUP_COLOURS

# heatmap colours: soft warm off-white = absent, deep wine = present
ABSENT, PRESENT = "#F5F1EC", "#7B1E3A"

meta, X = load()
n = len(meta)
bp_order = bioproject_order(meta)
bp_col = bioproject_colours(bp_order)

# isolate order: source group, then BioProject, then similarity (average linkage on Jaccard) inside each block
order = []
for g in GROUP_ORDER:
    for bp in bp_order:
        ids = list(meta.index[(meta["group"] == g) & (meta["bioproject"] == bp)])
        if len(ids) > 2:
            ids = [ids[i] for i in leaves_list(linkage(pdist(X.loc[ids].values.astype(bool), "jaccard"), "average"))]
        order += ids
genes = list(X.sum().sort_values(ascending=False).index)  # most frequent element first
M = X.loc[order, genes].values

fig = plt.figure(figsize=(16, 9.5))
gs = fig.add_gridspec(1, 3, width_ratios=[0.015, 0.015, 1], wspace=0.01, left=0.06, right=0.84, top=0.90, bottom=0.17)
ax_g, ax_b, ax = (fig.add_subplot(s) for s in gs)

ax.imshow(M, aspect="auto", cmap=ListedColormap([ABSENT, PRESENT]), interpolation="nearest")
ax.set_xticks(range(len(genes)))
ax.set_xticklabels(genes, rotation=90, fontsize=6)
ax.set_yticks([])
ax.tick_params(length=0)
# thin white lines between cells so dense blocks do not merge into one solid wall
ax.set_xticks(np.arange(-0.5, len(genes), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=0.35)
ax.tick_params(which="minor", length=0)
for sp in ax.spines.values():
    sp.set_visible(False)

for a, key, cols, title in ((ax_g, "group", GROUP_COLOURS, "Source"), (ax_b, "bioproject", bp_col, "BioProject")):
    a.imshow([[matplotlib.colors.to_rgb(cols[meta.loc[i, key]])] for i in order], aspect="auto", interpolation="nearest")
    a.set_xticks([]); a.set_yticks([])
    a.set_xlabel(title, fontsize=7, rotation=90, labelpad=4)
    for sp in a.spines.values():
        sp.set_visible(False)

# group labels (text, so the orange is never the only cue) and boundaries
pos = 0
for g in GROUP_ORDER:
    k = int((meta["group"] == g).sum())
    ax_g.text(-1.2, pos + k / 2, f"{g}\n(n = {k})", rotation=90, ha="right", va="center", fontsize=8)
    if pos:
        ax.axhline(pos - 0.5, color="white", lw=2)
    pos += k

handles = [Patch(color=GROUP_COLOURS[g], label=f"{g} (n = {(meta['group'] == g).sum()})") for g in GROUP_ORDER]
handles += [Patch(color=bp_col[b], label=f"BioProject: {b} (n = {(meta['bioproject'] == b).sum()})") for b in bp_order]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.855, 0.90), fontsize=7, frameon=False)
fig.suptitle(f"AMR genes and mutations in {n} E. coli isolates from Nigeria (Nextflow Bactopia)", fontsize=11, x=0.06, ha="left")
fig.text(0.06, 0.925, f"Each row is one isolate and each column is one AMR gene or mutation; coloured = detected. "
         f"n = {n} isolates, {len(genes)} elements. Rows are grouped by source, then BioProject, then similarity.", fontsize=7.5, ha="left")

Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/fig1_heatmap.pdf")
fig.savefig("figures/fig1_heatmap.png", dpi=300)
print("wrote figures/fig1_heatmap.pdf and .png;", n, "isolates,", len(genes), "elements")
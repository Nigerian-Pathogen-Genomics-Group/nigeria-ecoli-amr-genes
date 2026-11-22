#!/usr/bin/env python3
"""Issue #13, figures 2 and 3: PCoA of AMR gene profiles in the 117 E. coli.

Run from the repo root:
    PYTHONPATH=scripts python3 scripts/fig2_pcoa.py                        # dots coloured by source group
    PYTHONPATH=scripts python3 scripts/fig2_pcoa.py --colour-by bioproject  # dots coloured by BioProject
    PYTHONPATH=scripts python3 scripts/fig2_pcoa.py --colour-by both        # dots = BioProject, grey ellipses = source group
Writes figures/fig2_pcoa_source.* or figures/fig3_pcoa_bioproject.* (pdf and png),
results/pcoa_coordinates_117.tsv and results/pcoa_eigenvalues_117.tsv.

Method: Jaccard distance between the 0/1 AMR-element profiles of two isolates (1 - shared / present in either;
two isolates with no elements are given distance 0), then classical PCoA (double-centred squared distances,
eigendecomposition). Percent variance = eigenvalue / sum of the positive eigenvalues.
Isolates with identical profiles sit on the same point, so some points overlap.
Packages: numpy, pandas, matplotlib (print versions with: python3 -c "import numpy,pandas,matplotlib as m;print(numpy.__version__,pandas.__version__,m.__version__)")
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from fig_common import load, bioproject_order, bioproject_colours, GROUP_ORDER, GROUP_COLOURS

ap = argparse.ArgumentParser()
ap.add_argument("--colour-by", choices=["group", "bioproject", "both"], default="group")
args = ap.parse_args()

meta, X = load()
n = len(meta)
B = X.values.astype(bool)

# Jaccard distance matrix
inter = (B[:, None, :] & B[None, :, :]).sum(2)
union = (B[:, None, :] | B[None, :, :]).sum(2)
D = np.where(union == 0, 0.0, 1.0 - inter / np.maximum(union, 1))

# classical PCoA
J = np.eye(n) - np.ones((n, n)) / n
G = -0.5 * J @ (D ** 2) @ J
w, V = np.linalg.eigh(G)
idx = np.argsort(w)[::-1]
w, V = w[idx], V[:, idx]
pct = 100 * w / w[w > 1e-10].sum()
coords = V[:, :2] * np.sqrt(np.maximum(w[:2], 0))
n_neg = int((w < -1e-10).sum())
neg_share = 100 * (-w[w < -1e-10]).sum() / np.abs(w).sum()
print(f"isolates: {n}; elements: {B.shape[1]}")
print(f"PC1 {pct[0]:.1f}%  PC2 {pct[1]:.1f}%  (of the sum of positive eigenvalues)")
print(f"negative eigenvalues: {n_neg}; their share of the total absolute eigenvalue = {neg_share:.1f}%")

Path("results").mkdir(exist_ok=True)
out = meta[["group", "bioproject"]].copy()
out["PC1"], out["PC2"] = coords[:, 0], coords[:, 1]
out.to_csv("results/pcoa_coordinates_117.tsv", sep="\t")
pd.DataFrame({"axis": range(1, n + 1), "eigenvalue": w, "percent_of_positive_sum": pct}).to_csv(
    "results/pcoa_eigenvalues_117.tsv", sep="\t", index=False)

if args.colour_by == "group":
    series, keys, cols = meta["group"], GROUP_ORDER, GROUP_COLOURS
    name = lambda k: k
    stem = "fig2_pcoa_source"
    legend_kw = dict(loc="upper right")
else:
    series, keys = meta["bioproject"], bioproject_order(meta)
    cols = bioproject_colours(keys)
    name = (lambda k: k) if args.colour_by == "both" else (lambda k: f"BioProject: {k}")
    stem = "fig3_pcoa_bioproject" if args.colour_by == "bioproject" else "fig2_pcoa_both"
    legend_kw = dict(loc="upper left", bbox_to_anchor=(1.01, 1.0))
width = {"group": 8.5, "bioproject": 10.5, "both": 11.5}[args.colour_by]

fig, ax = plt.subplots(figsize=(width, 6.5))
ELL = {"Human": ("-", None), "Animal/food-chain": ((0, (6, 3)), None), "Environment": ((0, (1, 1.6)), None)}
if args.colour_by == "both":
    # one 68% data ellipse per source group (grey, different line style), drawn behind the BioProject dots
    for g, (ls, _) in ELL.items():
        xy = coords[(meta["group"] == g).values]
        mu, cov = xy.mean(0), np.cov(xy.T)
        t = np.linspace(0, 2 * np.pi, 361)
        pts = mu[:, None] + np.sqrt(2.279) * np.linalg.cholesky(cov) @ np.vstack([np.cos(t), np.sin(t)])
        ax.fill(pts[0], pts[1], facecolor="#888888", alpha=0.10, zorder=1, linewidth=0)
        ax.plot(pts[0], pts[1], color="#333333", lw=2.0, ls=ls, zorder=2)

for k in keys:
    m = (series == k).values
    ax.scatter(coords[m, 0], coords[m, 1], s=60, c=cols[k], alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
dots = [Line2D([], [], marker="o", ls="", markersize=8, markerfacecolor=cols[k], markeredgecolor="white",
               label=f"{name(k)} (n = {(series == k).sum()})") for k in keys]
if args.colour_by == "both":
    leg1 = ax.legend(handles=dots, title="BioProject", title_fontsize=10, frameon=False, fontsize=9.5,
                     alignment="left", loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.add_artist(leg1)
    lines = [Line2D([], [], color="#333333", lw=2.0, ls=ELL[g][0], label=f"{g} (n = {(meta['group'] == g).sum()})") for g in GROUP_ORDER]
    ax.legend(handles=lines, title="Source group", title_fontsize=10, frameon=False, fontsize=9.5,
              alignment="left", loc="upper left", bbox_to_anchor=(1.01, 0.50))
else:
    ax.legend(handles=dots, frameon=False, fontsize=9.5, **legend_kw)
ax.set_xlabel(f"PCoA 1 ({pct[0]:.1f}% of variance)")
ax.set_ylabel(f"PCoA 2 ({pct[1]:.1f}% of variance)")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.set_title(f"Resistance gene profiles of {n} E. coli isolates", fontsize=11, loc="left")
fig.text(0.01, 0.012, f"PCoA of Jaccard distances; n = {n} isolates, {B.shape[1]} AMR elements. "
         "Isolates with identical profiles overlap." + ("\nEllipses: each encloses about 68% of the isolates of one source group." if args.colour_by == "both" else ""), fontsize=7.5)
fig.tight_layout(rect=(0, 0.06 if args.colour_by == "both" else 0.035, 1, 1))
Path("figures").mkdir(exist_ok=True)
fig.savefig(f"figures/{stem}.pdf", bbox_inches="tight")
fig.savefig(f"figures/{stem}.png", dpi=300, bbox_inches="tight")
print(f"wrote figures/{stem}.pdf and .png")
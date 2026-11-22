#!/usr/bin/env python3
"""Issue #13, figure 3: how often each AMR element occurs in each source group (117 E. coli).

Run from the repo root:
    PYTHONPATH=scripts python3 scripts/fig3_frequency.py
Writes figures/fig3_frequency.pdf and .png and results/gene_frequency_by_source_117.tsv (all 86 elements).

Figure rule (fixed before looking at the result): show the elements found in at least MIN_ISOLATES of the 117 isolates;
the table has every element. Each dot = percent of the isolates of that source group that carry the element.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig_common import load, GROUP_ORDER, GROUP_COLOURS

MIN_ISOLATES = 5

meta, X = load()
n = len(meta)
size = meta["group"].value_counts().reindex(GROUP_ORDER)
count = pd.DataFrame({g: X[(meta["group"] == g).values].sum() for g in GROUP_ORDER})
pct = count / size * 100
tab = pd.DataFrame({"element": X.columns, "isolates_total": X.sum().values})
for g in GROUP_ORDER:
    tab[f"n_{g}"] = count[g].values
    tab[f"percent_{g}"] = pct[g].round(1).values
tab["percent_total"] = (X.sum().values / n * 100).round(1)
tab = tab.sort_values(["isolates_total", "element"], ascending=[False, True])
Path("results").mkdir(exist_ok=True)
tab.to_csv("results/gene_frequency_by_source_117.tsv", sep="\t", index=False)

show = tab[tab["isolates_total"] >= MIN_ISOLATES].reset_index(drop=True)
print(f"{n} isolates; {len(tab)} elements; {len(show)} shown (found in >= {MIN_ISOLATES} isolates)")

fig, ax = plt.subplots(figsize=(10, 0.3 * len(show) + 1.6))
y = np.arange(len(show))[::-1]
for yi in y:
    ax.axhline(yi, color="#E6E6E6", lw=0.6, zorder=0)
for g, off in zip(GROUP_ORDER, (0.24, 0.0, -0.24)):  # small vertical offset so dots at the same value stay visible
    ax.scatter(show[f"percent_{g}"], y + off, s=42, c=GROUP_COLOURS[g], edgecolor="white", linewidth=0.5, alpha=0.9, zorder=3)
ax.set_yticks(y)
ax.set_yticklabels([f"{e}  ({k})" for e, k in zip(show["element"], show["isolates_total"])], fontsize=8)
ax.set_xlim(-2, 102)
ax.set_ylim(-0.7, len(show) - 0.3)
ax.set_xlabel("Isolates carrying the element (%)")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(handles=[Line2D([], [], marker="o", ls="", markersize=8, markerfacecolor=GROUP_COLOURS[g], markeredgecolor="white",
                          label=f"{g} (n = {size[g]})") for g in GROUP_ORDER],
          title="Source group", title_fontsize=10, frameon=False, fontsize=9.5, alignment="left",
          loc="upper left", bbox_to_anchor=(1.02, 1.0))
ax.set_title(f"AMR genes and mutations by source group in {n} E. coli isolates", fontsize=11, loc="left")
fig.text(0.01, 0.005, f"Showing the {len(show)} of {len(tab)} resistance genes and mutations found in at least {MIN_ISOLATES} of the {n} isolates. "
         "The number in brackets is how many isolates carry it.\nEach dot is the percentage of isolates in that source group that carry the element.", fontsize=7.5)
fig.tight_layout(rect=(0, 0.035, 1, 1))
Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/fig3_frequency.pdf", bbox_inches="tight")
fig.savefig("figures/fig3_frequency.png", dpi=300, bbox_inches="tight")
print("wrote figures/fig3_frequency.pdf and .png and results/gene_frequency_by_source_117.tsv")
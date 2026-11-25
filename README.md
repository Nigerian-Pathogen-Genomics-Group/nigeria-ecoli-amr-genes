# Genomic Survey of AMR Genes in Nigerian *E. coli*

A project of the **Nigerian Pathogen Genomics Group (NPGG)** and the **Nigerian Nextflow community**.

Our goal is a standard, open, reproducible pipeline that any African lab can run on raw genomes to screen for antimicrobial resistance (AMR) genes and sequence types. We test it here on publicly available Nigerian *Escherichia coli* genomes.

## Aim

Describe the AMR gene profiles of public Nigerian *E. coli* genomes from Lagos, Ibadan and Edo State, across human clinical, animal/food-chain and environmental sources (including e-waste and abattoir samples).

Places, sources and submitting studies overlap, so we describe what is found and **do not claim differences between sources or places**.

## What we do

1. **Reproduce and standardise (Nextflow).** We run [Bactopia](https://bactopia.github.io/) (Nextflow) on the assembled genomes to confirm the species (*E. coli*, not *Shigella*), assign sequence types (MLST), and re-call AMR genes with one tool version (AMRFinderPlus).
2. **Analyse (Python or R).** Isolates-by-genes table, gene frequency by source, heatmap, ordination (PCoA), and a sensitivity check without the largest submitting study.
3. **Predict antimicrobial peptides (in progress).** [Macrel](https://github.com/BigDataBiology/macrel), a machine-learning tool, is run through nf-core/funcscan. Results are preliminary and will be added when ready. Progress is tracked in the Issues tab.

## Data

Public isolates from NCBI Pathogen Detection. See `data/README.md` for the source query, the three files (raw 888, cleaned 166, analysis set 125) and exactly how each was made. Never edit the raw file.

All data are public. We cite the original studies and BioProjects, and every number in the manuscript must be traceable to the files in `data/` or to a cited source.

## Repository layout

```
docs/      proposal, literature review, manuscript drafts
data/      raw, cleaned and analysis files
scripts/   analysis code and Nextflow commands
results/   tables and figures
```

## Contribute

We welcome collaborators. You can do everything in your web browser.

1. **Claim a task.** Open the Issues tab, choose a task and comment "I'll take this".
2. **Add your work.** Open the right folder (`docs/`, `scripts/` or `results/`), click **Add file**, then **Upload files** or **Create new file**. To change a file, open it and click the pencil icon.
3. **Commit to a new branch.** Choose **Create a new branch for this commit and start a pull request**, not "Commit directly to the main branch". Name the branch with your name and task, for example `amina-literature-review`.
4. **Open the pull request.** Add a short title and write "Closes #number" (the issue's number) in the description.
5. **Review.** Nkiruka reviews and merges. Please do not merge your own pull request, and do not commit large raw outputs (ask in the issue first).

Your contribution is counted when your work is merged and the email on your commit is linked to your GitHub account (Settings, then Emails). If you get stuck, comment on your issue.

## Contributors and authorship

Authors are listed by contribution to the work (design, analysis, writing, review), as agreed in the proposal. After your first merged contribution, Nkiruka adds your name and role. If your name is missing, comment on your pull request or message her.

| Name | Role |
|---|---|
| Nkiruka Cynthia Efenji | Coordinator |
| Peace Daniel | Literature Review |

## Status

Draft. Target: manuscript ready December 2026 / January 2027.

# Genomic Survey of AMR Genes in Nigerian E. coli

Part of the Nigerian Pathogen Genomics Group (NPGG) and the Nigerian Nextflow community.

## Aim
Describe the antimicrobial resistance (AMR) gene profiles of publicly available
Nigerian *Escherichia coli* genomes from Lagos, Ibadan and Edo State, and what is found
in each setting (human clinical, e-waste and abattoir). Because places, sources and
submitting studies overlap, we do not claim differences between sources or places.

This is Project 1 (Stage 1). Predicted antimicrobial peptides (Macrel) are Project 2,
in a separate repository, and start only if this project works.

## Data
Public isolates from NCBI Pathogen Detection. See `data/README.md` for the source query,
the three files (raw 888, cleaned 166, analysis set 125) and exactly how each was made.
Never edit the raw file.

## Workflow
- **Nextflow (Bactopia):** we run Bactopia on the assembled genomes to confirm the species
  (*E. coli*, not *Shigella*), assign sequence types (MLST), and re-call resistance genes
  with one tool version.
- **Analysis (Python or R):** isolates-by-genes table, gene frequency by source, heatmap, PCA,
  and a sensitivity check without the largest submitting study.
- All code is in `scripts/`, outputs in `results/`, and documents in `docs/`.

## Repository layout
- `docs/` proposal, literature review, manuscript drafts
- `data/` raw, cleaned and analysis files
- `scripts/` analysis code and Nextflow commands
- `results/` tables and figures

## How we work
1. Pick a task from the **Issues** tab and comment "I'll take this" to claim it.
2. Do the work on your own branch (see the guide below).
3. Open a **pull request** and link the issue (write "Closes #number").
4. Nkiruka reviews and merges. Please do not merge your own pull request.
5. Never commit large raw outputs. Ask in the issue first.

## How to submit a pull request (no GitHub experience needed)
You can do everything in your web browser.

1. **Claim a task.** Open the Issues tab, choose a task, and comment "I'll take this".
2. **Add your work.** On the repo's main page, open the right folder (`docs/`, `scripts/`
   or `results/`). Click **Add file**, then **Upload files** (to add a document) or
   **Create new file** (to write text). To change an existing file, open it and click
   the pencil icon.
3. **Commit to a new branch.** Click **Commit changes**. In the box that appears, choose
   **Create a new branch for this commit and start a pull request**. Do NOT choose
   "Commit directly to the main branch". Name the branch with your name and the task,
   for example `amina-literature-review`. Click **Propose changes**.
4. **Open the pull request.** Add a short title (for example "Literature review draft"),
   write "Closes #number" in the description (use the issue's number), and click
   **Create pull request**.
5. **Wait for review.** Nkiruka will comment or merge. If she asks for changes, open
   your pull request, click the **Files changed** tab, edit the file again on the same
   branch, and commit. The pull request updates itself.
6. **Check your name.** Your contribution is counted when your work is merged and the
   email on your commit is linked to your GitHub account (Settings, then Emails).

If you get stuck, comment on your issue or ask in the group.

## Contributors and authorship
Authors are listed by contribution to the work (design, analysis, writing, review),
as agreed in the proposal. Roles are agreed in the group before writing starts.

**Nkiruka maintains the list below** After your first merged
contribution, Nkiruka adds your name and role. If you think your name is missing,
comment on your pull request or message her.

| Name | Role |
|---|---|
| Nkiruka Cynthia Efenji | Coordinator |

## Data use
All data are public. We cite the original studies and BioProjects, and contact the
submitting authors where appropriate. Every number in the manuscript must be traceable
to the files in `data/` or to a cited source.

## Status
Draft. Target: manuscript ready December 2026 / January 2027.

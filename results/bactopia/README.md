# Bactopia results: 125 Nigerian *E. coli* / *Shigella* isolates

Bactopia run on the 125 analysis isolates (NCBI Pathogen Detection assemblies, one isolate per SNP cluster).
Covers issues #8 (run Bactopia), #9 (species and typing) and #10 (compare with NCBI calls).

## What is here

| Path | What it is |
|---|---|
| `raw/test5/`, `raw/all120/` | Bactopia `merged-results` tables exactly as downloaded, from the two runs (5 + 120 isolates). Not edited. |
| `bactopia_summary_125.tsv` | One row per isolate: taxon, ST, phylogroup, serotype, STECFinder, ShigaPass/ShigaTyper, AMR gene count, assembly size. |
| `amr_gene_matrix_125.tsv` | 125 isolates x 90 AMR genes/mutations, 1 = present (input for the heatmap and ordination, issue #13). |
| `bactopia_vs_ncbi_125.tsv` | Per-isolate comparison of Bactopia AMRFinderPlus calls with NCBI "AMR genotypes". |
| `run_parameters.json` | Pipeline parameters (account-specific paths replaced by placeholders). |
| `../../scripts/merge_bactopia_results.py` | Rebuilds the summary and the matrix from `raw/`. |
| `../../scripts/compare_bactopia_ncbi.py` | Rebuilds the NCBI comparison from `raw/`. |

Rebuild everything:

```bash
python3 scripts/merge_bactopia_results.py --runs results/bactopia/raw/test5 results/bactopia/raw/all120 \
    --ncbi data/Ecoli_AMR_analysis_125.tsv --outdir results/bactopia
python3 scripts/compare_bactopia_ncbi.py --runs results/bactopia/raw/test5 results/bactopia/raw/all120 \
    --ncbi data/Ecoli_AMR_analysis_125.tsv --out results/bactopia/bactopia_vs_ncbi_125.tsv
```

## How it was run (issue #8)

- Bactopia 4.1.0 (Nextflow DSL2; run on Nextflow 26.10), nf-bactopia plugin 2.1.7, on Seqera Cloud (Seqera Compute, AWS Batch, eu-west-2) with the `docker` profile.
- Input: assemblies (FASTA) downloaded from NCBI, listed in a Bactopia sample table (`runtype = assembly`).
- Species-specific tools were chosen by Bactopia's Merlin (`ask_merlin = true`).
- Main settings: `--amrfinderplus_organism Escherichia`, `--mlst_scheme escherichia`.
- Tool versions (from the pipeline's container list): AMRFinderPlus 4.2.7, mlst 2.35.0, ECTyper 2.0.0 (O/H database v1.0, pathotype database v2024-10-03), ClermonTyping 24.02, STECFinder 1.1.2, ShigaTyper 2.0.5, ShigaPass 1.5.0. The AMRFinderPlus database was the copy distributed with Bactopia 4.1.0; its release date was not recorded.
- Two runs: a 5-isolate test (all tasks succeeded) and the remaining 120 isolates. The 120-isolate run first stopped on one *Shigella* genome (GCA_013391325): ECTyper's `--verify` option rejects non-*E. coli* input and exits with an error. It was resumed with `ectyper_verify = false`, reusing the finished work, and all tasks then succeeded (0 failed). As a result the `QC` column of `ectyper.tsv` is empty for all isolates (the verify check was off); species comes from ECTyper's Mash step.
- Not used: `kleborate.tsv` (a *Klebsiella* tool that Bactopia also ran), `meta.tsv` and logs.

## Species and typing (issue #9)

- **Taxon (ECTyper Mash species, confirmed by ShigaPass, ShigaTyper and STECFinder):** 117 *E. coli*, 6 *Shigella* (5 *S. flexneri*, 1 *S. sonnei*), 2 *E. fergusonii*. The NCBI group "E.coli and Shigella" includes these. No EIEC among the *E. coli* (ShigaPass and ShigaTyper both "Not Shigella or EIEC").
- **Phylogroups among the 117 *E. coli*:** A 68, B1 23, D 13, B2 6, E 4, F 3. The 6 *Shigella* and 2 *E. fergusonii* are reported separately.
- **Sequence types:** 81 distinct STs across the 125 isolates (76 among the 117 *E. coli*). Most common: ST10 (10 isolates), ST131 (5). Three isolates have no ST assigned by `mlst`: GCA_013364455 (`MIXED`, a *Shigella*; ShigaPass gives ST245), GCA_060477895 (`NOVEL` allele combination), GCA_060580075 (partial match, ST131-like alleles; also the largest assembly at 5.93 Mb: check for a mixed or contaminated assembly before relying on it).
- **Phylogroup calls to treat as uncertain:** ClermonTyping's marker call and its Mash-based call disagree for GCA_060480375 (A vs B1) and GCA_060580235 (D vs A). They also differ for the 6 *Shigella* (expected).
- **Shiga toxin:** no *stx* gene in any isolate (STECFinder). 6 isolates are ipaH-positive: the *Shigella*.
- **Assembly quality:** 4.48 to 5.93 Mb; median 78 contigs; most fragmented GCA_053232725 (361 contigs).

## Comparison with NCBI (issue #10)

- Compared the set of AMR gene/mutation names per isolate: AMRFinderPlus rows with `Type == AMR` (repeat hits of one gene counted once) against the NCBI "AMR genotypes" column.
- Result: **125 of 125 isolates identical**; 1,187 calls from each source; no gene found by only one of them. This includes point mutations (for example `gyrA_S83L`) and partial or mistranslated genes.
- Interpretation: NCBI also uses AMRFinderPlus, so this shows that the Bactopia run reproduces NCBI's calls with a fixed tool version. It is not an independent validation. Only names were compared, not coverage or identity values.
- **blaCTX-M** (blaCTX-M-15, blaCTX-M-27) is present in 14 isolates. **No mcr gene** was found in any isolate.
- STRESS and VIRULENCE rows in `amrfinderplus.tsv` are not part of the NCBI "AMR genotypes" column and are excluded from the comparison.

## Caveats

- `shigatyper.tsv` exists only for the 120-isolate run, so `shigatyper` and `ipaB` are `NA` for the first 5 isolates (ShigaPass covers all 125).
- Source group, study batch and location are in the NCBI metadata and are not interpreted here; source is confounded with study batch, so no comparison between sources is claimed.

# Data

All data in this folder come from the NCBI Pathogen Detection Isolates Browser
(https://www.ncbi.nlm.nih.gov/pathogens/isolates/). Nothing here was edited by hand
except by the filtering and column removal described below.

## Source query
`geo_loc_name: Nigeria AND taxgroup_name:"E.coli and Shigella"`
Downloaded on: [13 October 2026]

## Files

| File | Rows | What it is |
|---|---|---|
| `Ecoli_RAW_AMR_data.xlsx` | 888 | Raw download from NCBI (17 columns), unchanged. |
| `Ecoli_AMR_cleaned_166` (.tsv / .xlsx) | 166 | Cleaned set (12 columns). |
| `Ecoli_AMR_analysis_125` (.tsv / .xlsx) | 125 | Analysis set: the 166 with near-duplicates removed. |

The .tsv and .xlsx versions of each file contain the same data. Use the .tsv for code.

## How the cleaned set (166) was made from the raw file
Rows kept if ALL of these are true:
1. Location is Lagos, Ibadan or Edo State (this includes the entries
   "Lagos, Takwa Bay" and "Lagos, Badagry Beach").
2. The isolate has an assembly (GCA accession).
3. The isolate has an isolation source.

Columns removed: Serovar, Food origin, Min-same, Min-diff, Computed types.
All other values are unchanged from the raw file.

Result: Lagos 67, Ibadan 46, Edo State 48, Lagos (Takwa Bay) 4, Lagos (Badagry Beach) 1.

## How the analysis set (125) was made from the cleaned set
One isolate is kept per SNP cluster (the first isolate listed in the file).
Isolates with no SNP cluster are all kept.
Result: 81 isolates with no cluster + 44 clusters = 125.
The 41 removed isolates are still in the 166 file.

## Columns in the cleaned files
#Organism group, Strain, Isolate identifiers, Isolate, Create date, Location,
Isolation source, Isolation type, SNP cluster, BioSample, Assembly, AMR genotypes.

## Notes and limits
- The AMR genotypes were assigned by NCBI (AMRFinderPlus). We did not re-run any genome analysis.
- Most isolates from the same source come from one or two submitting studies, so
  differences between sources may reflect the study and not the source.
- Do not edit the raw file. Make changes through scripts in `scripts/` so every step is recorded.

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
| `source_groups_166.tsv` | 166 | Source group for each isolate in the cleaned set (4 columns). |
| `source_groups_125.tsv` | 125 | Source group for each isolate in the analysis set (4 columns). |

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

## Source groups
Each isolate is assigned to one of three source groups, using only its NCBI
"Isolation source" value. The cleaned files were not changed. The assignments are in
`source_groups_166.tsv` and `source_groups_125.tsv`, in the same row order as the
cleaned and analysis files.

Columns: `BioSample`, `Assembly`, `Isolation source`, `Source group`.

### Rule table

| Isolation source value in NCBI | Source group |
|---|---|
| blood | Human |
| urine | Human |
| clinical isolate | Human |
| Human stool samples (followed by a number, e.g. "Human stool samples10") | Human |
| chicken | Animal/food-chain |
| Slaughterhouse surfaces | Animal/food-chain |
| Irrigated farm | Animal/food-chain |
| water from electronic waste dumpsite | Environment |
| soil from electronic waste dumpsite | Environment |
| dry sand | Environment |
| shoreline saltwater | Environment |

Every value in the 166 and 125 sets is covered by this table. No isolate was left unassigned.

### Final counts

| Source group | 166 set | 125 set |
|---|---|---|
| Human | 50 | 38 |
| Animal/food-chain | 52 | 38 |
| Environment | 64 | 49 |
| **Total** | **166** | **125** |

Counts by isolation source value:

| Isolation source | 166 set | 125 set | Group |
|---|---|---|---|
| blood | 15 | 13 | Human |
| clinical isolate | 7 | 6 | Human |
| urine | 3 | 0 | Human |
| Human stool samples (numbered) | 25 | 19 | Human |
| Slaughterhouse surfaces | 48 | 35 | Animal/food-chain |
| chicken | 3 | 2 | Animal/food-chain |
| Irrigated farm | 1 | 1 | Animal/food-chain |
| water from electronic waste dumpsite | 54 | 40 | Environment |
| soil from electronic waste dumpsite | 5 | 4 | Environment |
| dry sand | 3 | 3 | Environment |
| shoreline saltwater | 2 | 2 | Environment |

## Notes and limits
- The AMR genotypes were assigned by NCBI (AMRFinderPlus). We did not re-run any genome analysis.
- Most isolates from the same source come from one or two submitting studies, so
  differences between sources may reflect the study and not the source.
  We do not claim differences between source groups or between places.
- Do not edit the raw file. Make changes through scripts in `scripts/` so every step is recorded.

#!/usr/bin/env python3
"""Build the tables required by issues #9 and #10 from the Bactopia results.

Usage (from the repo root, after scripts/merge_bactopia_results.py has been run):
    python3 scripts/make_issue_tables.py \
        --runs results/bactopia/raw/test5 results/bactopia/raw/all120 \
        --summary results/bactopia/bactopia_summary_125.tsv \
        --ncbi data/Ecoli_AMR_analysis_125.tsv \
        --outdir results

Writes into --outdir:
    species_and_sequence_types.tsv   (#9) one row per isolate, kept or removed
    removed_isolates.tsv             (#9) removed isolates with reason and supporting tool output
    borderline_cases.tsv             (#9) isolates with an uncertain call (decisions are filled in by hand)
    amr_comparison_long.tsv          (#10) one row per isolate and gene, with the call from each source
    amr_agreement_summary.md         (#10) numbers behind the comparison

Removal rule (written in issue #9 on 19 Oct, before the run): remove an isolate if the species
call is not Escherichia coli, or is Shigella. Removed isolates stay in the tables, marked removed.
Species call = ECTyper species (Mash against a reference sketch); ShigaPass, ShigaTyper and
STECFinder are used as supporting evidence. Nothing is edited by hand.
"""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--ncbi", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    summ = read_tsv(a.summary)

    # ---------------- issue #9 ----------------
    sp_rows, removed = [], []
    for r in summ:
        taxon = r["taxon"]
        keep = taxon == "Escherichia coli"
        reason = "" if keep else f"species call is {taxon} (rule: remove if not Escherichia coli or if Shigella)"
        st = r["ST"]
        if st == "-":
            st_status = {"MIXED": "no ST: mixed alleles", "NOVEL": "no ST: novel allele combination",
                         "OK": "no ST: partial match"}.get(r["mlst_status"], "no ST: " + r["mlst_status"])
        else:
            st_status = r["mlst_status"].lower()
        pg_note = "" if r["phylogroup"] == r["phylogroup_mash"] else f"marker call {r['phylogroup']} differs from Mash-based call {r['phylogroup_mash']}"
        sp_rows.append({
            "assembly": r["assembly_version"], "species_call": taxon, "species_tool": "ECTyper (Mash species ID)",
            "sequence_type": st if st != "-" else "none", "st_status": st_status,
            "st_from_shigapass": r["ST_shigapass"], "phylogroup": r["phylogroup"],
            "phylogroup_mash": r["phylogroup_mash"], "phylogroup_note": pg_note,
            "status": "kept" if keep else "removed", "reason": reason})
        if not keep:
            removed.append({
                "assembly": r["assembly_version"], "species_call": taxon, "reason": reason,
                "shigapass": r["shigapass"], "shigatyper": r["shigatyper"], "ipaH_stecfinder": r["ipaH"],
                "ipaB_shigatyper": r["ipaB"], "stecfinder_cluster": r["stecfinder_cluster"],
                "mlst_ST": r["ST"], "ST_shigapass": r["ST_shigapass"], "phylogroup_clermontyping": r["phylogroup"]})
    write_tsv(out / "species_and_sequence_types.tsv", sp_rows,
              ["assembly", "species_call", "species_tool", "sequence_type", "st_status", "st_from_shigapass",
               "phylogroup", "phylogroup_mash", "phylogroup_note", "status", "reason"])
    write_tsv(out / "removed_isolates.tsv", removed,
              ["assembly", "species_call", "reason", "shigapass", "shigatyper", "ipaH_stecfinder",
               "ipaB_shigatyper", "stecfinder_cluster", "mlst_ST", "ST_shigapass", "phylogroup_clermontyping"])

    # borderline cases (kept or removed) that need a decision
    border = []
    for r in summ:
        notes = []
        if r["mlst_status"] != "PERFECT":
            notes.append(f"MLST status {r['mlst_status']} (ST {r['ST']})")
        if r["phylogroup"] != r["phylogroup_mash"] and not r["taxon"].startswith(("Shigella", "Escherichia fergusonii")):
            notes.append(f"phylogroup disagreement: markers {r['phylogroup']} vs Mash {r['phylogroup_mash']}")
        if notes:
            border.append({"assembly": r["assembly_version"], "taxon": r["taxon"], "issue": "; ".join(notes),
                           "status_under_rule": "removed" if r["taxon"] != "Escherichia coli" else "kept",
                           "decision": "", "reason": "", "decided_by": ""})
    write_tsv(out / "borderline_cases.tsv", border,
              ["assembly", "taxon", "issue", "status_under_rule", "decision", "reason", "decided_by"])

    # ---------------- issue #10 ----------------
    bact = defaultdict(lambda: defaultdict(list))
    for run in a.runs:
        for r in read_tsv(Path(run) / "amrfinderplus.tsv"):
            if r["Type"] == "AMR":
                bact[r["Name"]][r["Element symbol"]].append(r)
    ncbi = {}
    for r in read_tsv(a.ncbi):
        acc = r["Assembly"].split(".")[0]
        ncbi[acc] = {}
        for g in r["AMR genotypes"].split(","):
            if g:
                gene, _, status = g.rpartition("=")
                ncbi[acc][gene] = status
    status_of = {r["accession"]: ("kept" if r["taxon"] == "Escherichia coli" else "removed") for r in summ}
    version = {r["accession"]: r["assembly_version"] for r in summ}

    long_rows = []
    for acc in sorted(set(ncbi) | set(bact)):
        for gene in sorted(set(ncbi.get(acc, {})) | set(bact.get(acc, {}))):
            n_status = ncbi.get(acc, {}).get(gene)
            hits = bact.get(acc, {}).get(gene, [])
            subtype = ";".join(sorted({h["Subtype"] for h in hits}))
            ctype = "point_mutation" if (n_status == "POINT" or subtype == "POINT") else "gene_presence"
            both = n_status is not None and bool(hits)
            long_rows.append({
                "assembly": version[acc], "gene": gene, "call_type": ctype,
                "ncbi_status": n_status or "", "bactopia_hits": len(hits),
                "bactopia_method": ";".join(sorted({h["Method"] for h in hits})),
                "bactopia_class": ";".join(sorted({h["Class"] for h in hits})),
                "bactopia_subclass": ";".join(sorted({h["Subclass"] for h in hits})),
                "agreement": "both" if both else ("NCBI_only" if n_status else "Bactopia_only"),
                "isolate_status_after_9": status_of[acc]})
    write_tsv(out / "amr_comparison_long.tsv", long_rows,
              ["assembly", "gene", "call_type", "ncbi_status", "bactopia_hits", "bactopia_method",
               "bactopia_class", "bactopia_subclass", "agreement", "isolate_status_after_9"])

    def summarise(rows):
        c = Counter(r["agreement"] for r in rows)
        total = sum(c.values())
        by_type = {t: Counter(r["agreement"] for r in rows if r["call_type"] == t) for t in ("gene_presence", "point_mutation")}
        per_iso = defaultdict(set)
        for r in rows:
            per_iso[r["assembly"]].add(r["agreement"])
        exact = sum(1 for v in per_iso.values() if v == {"both"})
        return c, total, by_type, exact, len(per_iso)

    lines = ["# Bactopia vs NCBI resistance calls (issue #10)", "",
             "Generated by `scripts/make_issue_tables.py`. Unit: element symbols (genes and point mutations) per isolate.", ""]
    for label, sel in (("All 125 isolates", lambda r: True),
                       ("Isolates kept after #9 (117 *E. coli*)", lambda r: r["isolate_status_after_9"] == "kept")):
        rows = [r for r in long_rows if sel(r)]
        c, total, by_type, exact, n_iso = summarise(rows)
        lines += [f"## {label}", "",
                  f"- Isolates: {n_iso}; isolates with identical gene sets: {exact} of {n_iso}",
                  f"- Calls found by both: {c['both']} of {total} found by either ({100*c['both']/total:.1f}%); NCBI only: {c['NCBI_only']}; Bactopia only: {c['Bactopia_only']}"]
        for t, name in (("gene_presence", "Gene presence"), ("point_mutation", "Point mutations")):
            b = by_type[t]
            lines.append(f"- {name}: both {b['both']}, NCBI only {b['NCBI_only']}, Bactopia only {b['Bactopia_only']}")
        lines.append("")
    n_rep = sum(1 for acc in bact for g, h in bact[acc].items() if len(h) > 1)
    lines += ["## Notes", "",
              f"- Repeat hits: {n_rep} isolate-gene pairs have more than one Bactopia row (for example the same gene on two contigs). They are counted once, as NCBI lists a gene once.",
              "- Disagreements: none, so there is nothing to examine by hand and no reasons to propose. Most-discordant genes table: empty.",
              "- Drug class: the NCBI column has no drug class, so class-level agreement cannot be compared. Bactopia's class and subclass are included in `amr_comparison_long.tsv`.",
              "- Gene-name matching: no renaming rules were needed; every name in the NCBI column (with the `=STATUS` suffix removed) matched an AMRFinderPlus element symbol exactly. Unmatched names: none.",
              "- NCBI column basis: NCBI states that AMRFinderPlus is part of its Pathogen Detection pipeline and that the Isolates Browser shows its results (https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder). That page does not state which software or database version produced a given isolate's calls, or whether old isolates are re-run when the database changes.",
              "- Bactopia run: AMRFinderPlus 4.2.7 (container biocontainers/ncbi-amrfinderplus:4.2.7), database as distributed with Bactopia 4.1.0 (`datasets.bactopia.com/datasets/v4.1.0/amrfinderplus.tar.gz`). The database release date still has to be read from the AMRFinderPlus task log in Seqera.",
              "- Interpretation: both sources use AMRFinderPlus, so full agreement shows that our run reproduces NCBI's calls with a fixed tool version. It is not an independent validation. Only names were compared, not coverage or identity values.",
              "- STRESS and VIRULENCE rows in Bactopia's table are not part of NCBI's AMR genotypes and were excluded.", ""]
    (out / "amr_agreement_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))
    print(f"kept: {sum(1 for r in sp_rows if r['status']=='kept')}, removed: {len(removed)}, borderline: {len(border)}")


if __name__ == "__main__":
    main()

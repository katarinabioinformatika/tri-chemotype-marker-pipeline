#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re
import sys

BASE_DIR = Path(__file__).resolve().parents[1]

BLAST_DIR = BASE_DIR / "results" / "tri_blast"

META_FILE = BASE_DIR / "metadata" / "metadata.csv"

OUT_DIR = BASE_DIR / "results" / "tri_hits"

REFERENCE_HITS_FILE = OUT_DIR / "reconstructed_reference_hits.csv"

GENE_HITS_FILE = OUT_DIR / "tri_gene_hits.csv"

LOW_COVERAGE_FILE = OUT_DIR / "low_coverage_reference_hits.csv"

MERGE_DISTANCE_BP = 500

REFERENCE_COLLAPSE_DISTANCE_BP = 500

MIN_COMBINED_COVERAGE = 0.70

BLAST_COLUMNS = [
    "Reference_ID",
    "Chromosome_or_Contig",
    "Percent_identity",
    "Alignment_length",
    "Query_length",
    "Query_start",
    "Query_end",
    "Subject_start",
    "Subject_end",
    "Evalue",
    "Bitscore",
]

def extract_gene(reference_id):
    """Extract a standardised gene name from a reference sequence identifier."""

    first = str(reference_id).split("|")[0]

    match = re.search(
        r"(TRI\d+|PKS\d+|ZEB\d+)",
        first,
        re.IGNORECASE
    )

    if match:
        return match.group(1).upper()

    return first.upper()

def genome_id_from_filename(path):
    """Derive the genome identifier from a BLAST output filename."""

    name = path.name.replace("_tri_hits.tsv", "")

    return name.replace("_genomic", "")

def interval_union_length(intervals):
    """Calculate the number of unique positions covered by a set of intervals."""

    if not intervals:
        return 0

    intervals = [
        (min(start, end), max(start, end))
        for start, end in intervals
    ]

    intervals = sorted(intervals)

    merged = []

    for start, end in intervals:

        if not merged:
            merged.append([start, end])
            continue

        previous_start, previous_end = merged[-1]

        if start <= previous_end + 1:
            merged[-1][1] = max(previous_end, end)

        else:
            merged.append([start, end])

    total = sum(
        end - start + 1
        for start, end in merged
    )

    return total

def assign_reference_hit_groups(df):
    """Assign nearby HSPs from one reference sequence to reconstructed hits."""

    df = df.sort_values(
        ["Start", "End"]
    ).copy()

    groups = []

    current_group = 0

    current_end = None

    for _, row in df.iterrows():

        if current_end is None:
            current_group += 1
            current_end = row["End"]

        elif row["Start"] <= current_end + MERGE_DISTANCE_BP:

            current_end = max(
                current_end,
                row["End"]
            )

        else:
            current_group += 1
            current_end = row["End"]

        groups.append(current_group)

    df["Reference_hit_group"] = groups

    return df

def reconstruct_reference_hit(df):
    """Summarise several HSPs as one reconstructed reference hit."""

    query_length = int(
        df["Query_length"].iloc[0]
    )

    query_intervals = list(
        zip(
            df["Query_start"],
            df["Query_end"]
        )
    )

    query_covered_bp = interval_union_length(
        query_intervals
    )

    combined_coverage = (
        query_covered_bp / query_length
        if query_length > 0
        else 0
    )

    total_alignment_length = df[
        "Alignment_length"
    ].sum()

    if total_alignment_length > 0:

        weighted_identity = (
            (
                df["Percent_identity"]
                * df["Alignment_length"]
            ).sum()
            / total_alignment_length
        )

    else:
        weighted_identity = 0

    total_bitscore = df[
        "Bitscore"
    ].sum()

    best_evalue = df[
        "Evalue"
    ].min()

    return {
        "Genome_ID": df["Genome_ID"].iloc[0],
        "Gene": df["Gene"].iloc[0],
        "Reference_ID": df["Reference_ID"].iloc[0],
        "Chromosome_or_Contig": df[
            "Chromosome_or_Contig"
        ].iloc[0],
        "Strand": df["Strand"].iloc[0],
        "Reference_hit_group": df[
            "Reference_hit_group"
        ].iloc[0],

        "Start": int(df["Start"].min()),
        "End": int(df["End"].max()),

        "HSP_count": len(df),

        "Query_length": query_length,
        "Query_covered_bp": query_covered_bp,
        "Combined_query_coverage": combined_coverage,

        "Weighted_percent_identity": weighted_identity,
        "Total_bitscore": total_bitscore,
        "Best_evalue": best_evalue,

        "Total_alignment_length": int(
            total_alignment_length
        ),
    }

def assign_gene_locus_groups(df):
    """Group reconstructed alternative-reference hits at the same gene locus."""

    df = df.sort_values(
        ["Start", "End"]
    ).copy()

    groups = []

    current_group = 0

    current_end = None

    for _, row in df.iterrows():

        if current_end is None:
            current_group += 1
            current_end = row["End"]

        elif (
            row["Start"]
            <= current_end + REFERENCE_COLLAPSE_DISTANCE_BP
        ):

            current_end = max(
                current_end,
                row["End"]
            )

        else:
            current_group += 1
            current_end = row["End"]

        groups.append(current_group)

    df["Gene_locus_group"] = groups

    return df

def summarise_gene_locus(df):
    """Collapse alternative reference hits into one genomic gene occurrence."""

    ranked = df.sort_values(
        by=[
            "Combined_query_coverage",
            "Total_bitscore",
            "Weighted_percent_identity",
            "Best_evalue",
        ],
        ascending=[
            False,
            False,
            False,
            True,
        ]
    )

    representative = ranked.iloc[0]

    supporting_references = sorted(
        df["Reference_ID"]
        .astype(str)
        .unique()
    )

    total_supporting_hsps = int(
        df["HSP_count"].sum()
    )

    return {
        "Genome_ID": representative["Genome_ID"],
        "Gene": representative["Gene"],
        "Chromosome_or_Contig":
            representative["Chromosome_or_Contig"],
        "Strand": representative["Strand"],
        "Gene_locus_group":
            representative["Gene_locus_group"],

        "Start": int(df["Start"].min()),
        "End": int(df["End"].max()),

        "Representative_reference":
            representative["Reference_ID"],

        "Supporting_references":
            ";".join(supporting_references),

        "Number_of_supporting_references":
            len(supporting_references),

        "Number_of_reference_hits":
            len(df),

        "Total_supporting_HSPs":
            total_supporting_hsps,

        "Representative_HSP_count":
            int(representative["HSP_count"]),

        "Query_length":
            int(representative["Query_length"]),
        "Query_covered_bp":
            int(representative["Query_covered_bp"]),
        "Combined_query_coverage":
            representative["Combined_query_coverage"],

        "Weighted_percent_identity":
            representative["Weighted_percent_identity"],
        "Total_bitscore":
            representative["Total_bitscore"],
        "Best_evalue":
            representative["Best_evalue"],
    }

def main():
    """Reconstruct candidate genomic gene hits from the raw BLAST HSPs."""

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    blast_files = sorted(
        BLAST_DIR.glob("*_tri_hits.tsv")
    )

    if not blast_files:
        sys.exit(
            f"ERROR: No BLAST files found in {BLAST_DIR}"
        )

    all_hits = []

    for file in blast_files:

        genome_id = genome_id_from_filename(
            file
        )

        if file.stat().st_size == 0:
            print(
                f"WARNING: Empty file skipped: {file.name}"
            )
            continue

        df = pd.read_csv(
            file,
            sep="\t",
            header=None,
            names=BLAST_COLUMNS
        )

        df["Genome_ID"] = genome_id

        all_hits.append(df)

    if not all_hits:
        sys.exit(
            "ERROR: No BLAST hits loaded."
        )

    hits = pd.concat(
        all_hits,
        ignore_index=True
    )

    raw_hsp_count = len(hits)

    hits["Gene"] = hits[
        "Reference_ID"
    ].apply(extract_gene)

    hits["Start"] = hits[
        ["Subject_start", "Subject_end"]
    ].min(axis=1)

    hits["End"] = hits[
        ["Subject_start", "Subject_end"]
    ].max(axis=1)

    hits["Strand"] = hits.apply(
        lambda row: "+"
        if row["Subject_start"] <= row["Subject_end"]
        else "-",
        axis=1
    )

    hits["HSP_query_coverage"] = (
        hits["Alignment_length"]
        / hits["Query_length"]
    )

    grouped_hsps = []

    for _, sub in hits.groupby(
        [
            "Genome_ID",
            "Gene",
            "Reference_ID",
            "Chromosome_or_Contig",
            "Strand",
        ],
        sort=False
    ):

        grouped_hsps.append(
            assign_reference_hit_groups(sub)
        )

    hits = pd.concat(
        grouped_hsps,
        ignore_index=True
    )

    reconstructed_records = []

    for _, sub in hits.groupby(
        [
            "Genome_ID",
            "Gene",
            "Reference_ID",
            "Chromosome_or_Contig",
            "Strand",
            "Reference_hit_group",
        ],
        sort=False
    ):

        reconstructed_records.append(
            reconstruct_reference_hit(sub)
        )

    reference_hits = pd.DataFrame(
        reconstructed_records
    )

    reference_hits[
        "Combined_query_coverage_percent"
    ] = (
        reference_hits["Combined_query_coverage"]
        * 100
    )

    reference_hits["Pass_coverage_filter"] = (
        reference_hits["Combined_query_coverage"]
        >= MIN_COMBINED_COVERAGE
    )

    reference_hits.to_csv(
        REFERENCE_HITS_FILE,
        index=False
    )

    low_coverage = reference_hits[
        ~reference_hits["Pass_coverage_filter"]
    ].copy()

    low_coverage.to_csv(
        LOW_COVERAGE_FILE,
        index=False
    )

    retained = reference_hits[
        reference_hits["Pass_coverage_filter"]
    ].copy()

    if retained.empty:

        print("Done.")
        print(
            f"Raw BLAST HSPs loaded: {raw_hsp_count}"
        )
        print(
            f"Reconstructed reference hits: {len(reference_hits)}"
        )
        print(
            "No reconstructed hits passed the "
            f"{MIN_COMBINED_COVERAGE:.0%} coverage threshold."
        )
        print(
            f"Reference hits saved to: {REFERENCE_HITS_FILE}"
        )
        print(
            f"Low-coverage hits saved to: {LOW_COVERAGE_FILE}"
        )

        return

    locus_grouped = []

    for _, sub in retained.groupby(
        [
            "Genome_ID",
            "Gene",
            "Chromosome_or_Contig",
            "Strand",
        ],
        sort=False
    ):

        locus_grouped.append(
            assign_gene_locus_groups(sub)
        )

    retained = pd.concat(
        locus_grouped,
        ignore_index=True
    )

    gene_records = []

    for _, sub in retained.groupby(
        [
            "Genome_ID",
            "Gene",
            "Chromosome_or_Contig",
            "Strand",
            "Gene_locus_group",
        ],
        sort=False
    ):

        gene_records.append(
            summarise_gene_locus(sub)
        )

    gene_hits = pd.DataFrame(
        gene_records
    )

    gene_hits = gene_hits.sort_values(
        by=[
            "Genome_ID",
            "Gene",
            "Chromosome_or_Contig",
            "Start",
            "End",
        ]
    ).reset_index(drop=True)

    gene_hits["Gene_hit_number"] = (
        gene_hits.groupby(
            ["Genome_ID", "Gene"]
        ).cumcount()
        + 1
    )

    gene_hits["Gene_hit_ID"] = (
        gene_hits["Gene"]
        + "_hit_"
        + gene_hits["Gene_hit_number"].astype(str)
    )

    gene_hits["Multiple_hits_for_gene"] = (
        gene_hits.groupby(
            ["Genome_ID", "Gene"]
        )["Gene_hit_ID"]
        .transform("count")
        > 1
    )

    gene_hits["Multiple_contigs_for_gene"] = (
        gene_hits.groupby(
            ["Genome_ID", "Gene"]
        )["Chromosome_or_Contig"]
        .transform("nunique")
        > 1
    )

    gene_hits[
        "Combined_query_coverage_percent"
    ] = (
        gene_hits["Combined_query_coverage"]
        * 100
    )

    if META_FILE.exists():

        meta = pd.read_csv(
            META_FILE
        )

        keep_cols = [
            column
            for column in [
                "Genome_ID",
                "Species",
                "Chemotype",
                "Toxin_System",
                "Strain",
            ]
            if column in meta.columns
        ]

        gene_hits = gene_hits.merge(
            meta[keep_cols].drop_duplicates(),
            on="Genome_ID",
            how="left"
        )

    else:

        gene_hits["Species"] = ""
        gene_hits["Chemotype"] = ""
        gene_hits["Toxin_System"] = ""
        gene_hits["Strain"] = ""

    final_cols = [
        "Genome_ID",
        "Species",
        "Chemotype",
        "Toxin_System",

        "Gene",

        "Chromosome_or_Contig",
        "Start",
        "End",
        "Strand",

        "Multiple_hits_for_gene",
        "Multiple_contigs_for_gene",

        "Representative_HSP_count",

        "Query_length",
        "Query_covered_bp",
        "Combined_query_coverage_percent",

        "Weighted_percent_identity",
        "Total_bitscore",
        "Best_evalue",
    ]

    gene_hits = gene_hits[
        [
            column
            for column in final_cols
            if column in gene_hits.columns
        ]
    ]

    gene_hits.to_csv(
        GENE_HITS_FILE,
        index=False
    )

    print("\nDone.")

    print(
        f"Raw BLAST HSPs loaded: {raw_hsp_count}"
    )

    print(
        "Reconstructed reference hits: "
        f"{len(reference_hits)}"
    )

    print(
        "Reference hits passing "
        f"{MIN_COMBINED_COVERAGE:.0%} coverage: "
        f"{len(retained)}"
    )

    print(
        "Reference hits below "
        f"{MIN_COMBINED_COVERAGE:.0%} coverage: "
        f"{len(low_coverage)}"
    )

    print(
        "Final genomic gene hits: "
        f"{len(gene_hits)}"
    )

    print(
        f"Reconstructed reference hits saved to: "
        f"{REFERENCE_HITS_FILE}"
    )

    print(
        f"Low-coverage hits saved to: "
        f"{LOW_COVERAGE_FILE}"
    )

    print(
        f"Final gene hits saved to: "
        f"{GENE_HITS_FILE}"
    )

if __name__ == "__main__":
    main()
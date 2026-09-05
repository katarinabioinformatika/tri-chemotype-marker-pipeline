#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys

BASE_DIR = Path(__file__).resolve().parents[1]

IN_FILE = BASE_DIR / "results" / "tri_hits" / "tri_gene_hits.csv"

OUT_DIR = BASE_DIR / "results" / "tri_loci"

LOCI_FILE = OUT_DIR / "tri_loci.csv"

MEMBERS_FILE = OUT_DIR / "tri_locus_members.csv"

# Maximum gap between neighbouring hits in a locus; selected from the
# observed separation of TRI-gene hits in this dataset.
LOCUS_DISTANCE_BP = 15_000

def assign_loci(df):
    """Assign spatial locus numbers to candidate gene hits on one genomic sequence."""

    df = df.sort_values(
        ["Start", "End"]
    ).copy()

    locus_numbers = []

    current_locus = 0

    current_locus_end = None

    for _, row in df.iterrows():

        start = int(row["Start"])
        end = int(row["End"])

        if current_locus_end is None:
            current_locus += 1
            current_locus_end = end

        else:
            gap_from_current_locus = (
                start - current_locus_end - 1
            )

            if gap_from_current_locus <= LOCUS_DISTANCE_BP:

                current_locus_end = max(
                    current_locus_end,
                    end
                )

            else:
                current_locus += 1
                current_locus_end = end

        locus_numbers.append(
            current_locus
        )

    df["Local_Locus_Number"] = locus_numbers

    df["Previous_End"] = (
        df.groupby(
            "Local_Locus_Number"
        )["End"]
        .shift(1)
    )

    df["Next_Start"] = (
        df.groupby(
            "Local_Locus_Number"
        )["Start"]
        .shift(-1)
    )

    df["Gap_to_previous"] = (
        df["Start"]
        - df["Previous_End"]
        - 1
    )

    df["Gap_to_next"] = (
        df["Next_Start"]
        - df["End"]
        - 1
    )

    df = df.drop(
        columns=[
            "Previous_End",
            "Next_Start",
        ]
    )

    return df

def summarise_loci(members):
    """Summarise candidate gene members as one record per spatial locus."""

    summary = (
        members
        .groupby(
            [
                "Genome_ID",
                "Chromosome_or_Contig",
                "Locus_ID",
            ],
            as_index=False
        )
        .agg(
            Locus_start=(
                "Start",
                "min"
            ),
            Locus_end=(
                "End",
                "max"
            ),

            Number_of_hits=(
                "Gene",
                "count"
            ),

            Number_of_unique_genes=(
                "Gene",
                "nunique"
            ),

            Genes=(
                "Gene",
                lambda x: ",".join(
                    x.astype(str)
                )
            ),

            Maximum_gap_bp=(
                "Gap_to_next",
                "max"
            ),
        )
    )

    summary["Locus_length_bp"] = (
        summary["Locus_end"]
        - summary["Locus_start"]
        + 1
    )

    summary["Maximum_gap_bp"] = (
        summary["Maximum_gap_bp"]
        .fillna(0)
        .astype(int)
    )

    return summary

def main():
    """Group retained candidate genomic gene hits into spatial loci."""

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    if not IN_FILE.exists():
        sys.exit(
            f"ERROR: Input file not found: {IN_FILE}"
        )

    hits = pd.read_csv(
        IN_FILE
    )

    if hits.empty:
        sys.exit(
            f"ERROR: Input file contains no gene hits: {IN_FILE}"
        )

    required_cols = [
        "Genome_ID",
        "Chromosome_or_Contig",
        "Gene",
        "Start",
        "End",
    ]

    missing = [
        column
        for column in required_cols
        if column not in hits.columns
    ]

    if missing:
        sys.exit(
            f"ERROR: Missing required columns: {missing}"
        )

    hits["Start"] = hits["Start"].astype(int)
    hits["End"] = hits["End"].astype(int)

    grouped = []

    for _, sub in hits.groupby(
        [
            "Genome_ID",
            "Chromosome_or_Contig",
        ],
        sort=False
    ):

        grouped.append(
            assign_loci(sub)
        )

    members = pd.concat(
        grouped,
        ignore_index=True
    )

    members["Locus_ID"] = (
        members["Genome_ID"].astype(str)
        + "|"
        + members["Chromosome_or_Contig"].astype(str)
        + "|L"
        + members["Local_Locus_Number"].astype(str)
    )

    members = members.sort_values(
        by=[
            "Genome_ID",
            "Chromosome_or_Contig",
            "Start",
            "End",
        ]
    ).reset_index(
        drop=True
    )

    locus_summary = summarise_loci(
        members
    )

    metadata_cols = [
        column
        for column in [
            "Genome_ID",
            "Species",
            "Chemotype",
            "Toxin_System",
            "Strain",
        ]
        if column in members.columns
    ]

    if len(metadata_cols) > 1:

        metadata = (
            members[metadata_cols]
            .drop_duplicates(
                subset=["Genome_ID"]
            )
        )

        locus_summary = locus_summary.merge(
            metadata,
            on="Genome_ID",
            how="left"
        )

    locus_cols = [
        "Genome_ID",
        "Species",
        "Chemotype",
        "Toxin_System",
        "Strain",

        "Chromosome_or_Contig",
        "Locus_ID",

        "Locus_start",
        "Locus_end",
        "Locus_length_bp",

        "Number_of_hits",
        "Number_of_unique_genes",
        "Maximum_gap_bp",

        "Genes",
    ]

    locus_summary = locus_summary[
        [
            column
            for column in locus_cols
            if column in locus_summary.columns
        ]
    ]

    locus_summary = locus_summary.sort_values(
        by=[
            "Genome_ID",
            "Chromosome_or_Contig",
            "Locus_start",
            "Locus_end",
        ]
    ).reset_index(
        drop=True
    )

    member_cols = [
        "Genome_ID",
        "Species",
        "Chemotype",
        "Toxin_System",
        "Strain",

        "Chromosome_or_Contig",
        "Locus_ID",
        "Local_Locus_Number",

        "Gene",
        "Start",
        "End",
        "Strand",

        "Gap_to_previous",
        "Gap_to_next",

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

    members = members[
        [
            column
            for column in member_cols
            if column in members.columns
        ]
    ]

    locus_summary.to_csv(
        LOCI_FILE,
        index=False
    )

    members.to_csv(
        MEMBERS_FILE,
        index=False
    )

    print("\nDone.")

    print(
        f"Input gene hits: {len(hits)}"
    )

    print(
        f"Spatial loci identified: {len(locus_summary)}"
    )

    print(
        "Maximum neighbouring-hit gap: "
        f"{LOCUS_DISTANCE_BP:,} bp"
    )

    print(
        f"Locus summary saved to: {LOCI_FILE}"
    )

    print(
        f"Locus members saved to: {MEMBERS_FILE}"
    )

if __name__ == "__main__":
    main()

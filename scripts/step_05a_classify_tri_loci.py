#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys

BASE_DIR = Path(__file__).resolve().parents[1]

IN_FILE = BASE_DIR / "results" / "tri_loci" / "tri_loci.csv"

OUT_DIR = BASE_DIR / "results" / "classified_loci"

OUT_FILE = OUT_DIR / "tri_loci_classified.csv"

# TRI5 anchors the core cluster, but two additional core genes are required
# to prevent an isolated TRI5 hit from being classified as a complete cluster.
CORE_TRI_GENES = {
    "TRI3",
    "TRI4",
    "TRI5",
    "TRI6",
    "TRI7",
    "TRI8",
    "TRI9",
    "TRI10",
    "TRI11",
    "TRI12",
    "TRI13",
    "TRI14",
}

MIN_CORE_TRI_GENES = 3

def parse_genes(gene_string):
    """Convert the comma-separated Step 4 Genes field into a standardised set."""

    if pd.isna(gene_string):
        return set()

    return {
        gene.strip().upper()
        for gene in str(gene_string).split(",")
        if gene.strip()
    }

def classify_locus(gene_string):
    """Assign a biological TRI-locus class according to gene composition."""

    genes = parse_genes(
        gene_string
    )

    core_count = len(
        genes & CORE_TRI_GENES
    )

    if (
        "TRI5" in genes
        and core_count >= MIN_CORE_TRI_GENES
    ):
        return "core_TRI_cluster"

    elif (
        "TRI1" in genes
        or "TRI16" in genes
    ):
        return "TRI1_TRI16_locus"

    elif "TRI101" in genes:
        return "TRI101_locus"

    else:
        return "other_TRI_locus"

def main():
    """Classify spatial TRI loci according to their gene composition."""

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    if not IN_FILE.exists():
        sys.exit(
            f"ERROR: Input file not found: {IN_FILE}"
        )

    loci = pd.read_csv(
        IN_FILE
    )

    if loci.empty:
        sys.exit(
            f"ERROR: Input file contains no loci: {IN_FILE}"
        )

    if "Genes" not in loci.columns:
        sys.exit(
            "ERROR: Input file must contain a 'Genes' column."
        )

    loci["Core_TRI_gene_count"] = (
        loci["Genes"]
        .apply(
            lambda gene_string:
            len(
                parse_genes(gene_string)
                & CORE_TRI_GENES
            )
        )
    )

    loci["TRI5_present"] = (
        loci["Genes"]
        .apply(
            lambda gene_string:
            "TRI5" in parse_genes(gene_string)
        )
    )

    loci["Locus_class"] = (
        loci["Genes"]
        .apply(
            classify_locus
        )
    )

    preferred_cols = [
        "Genome_ID",
        "Species",
        "Chemotype",
        "Toxin_System",
        "Strain",

        "Chromosome_or_Contig",
        "Locus_ID",

        "Locus_class",

        "Locus_start",
        "Locus_end",
        "Locus_length_bp",

        "Number_of_hits",
        "Number_of_unique_genes",

        "TRI5_present",
        "Core_TRI_gene_count",

        "Maximum_gap_bp",

        "Genes",
    ]

    loci = loci[
        [
            column
            for column in preferred_cols
            if column in loci.columns
        ]
    ]

    loci.to_csv(
        OUT_FILE,
        index=False
    )

    print("\nDone.")

    print(
        f"Input loci: {len(loci)}"
    )

    print(
        "Core TRI cluster criterion: "
        f"TRI5 present + at least "
        f"{MIN_CORE_TRI_GENES} core TRI genes"
    )

    print(
        f"Output saved to: {OUT_FILE}"
    )

    print("\nLocus class counts:")

    print(
        loci["Locus_class"]
        .value_counts()
        .to_string()
    )

if __name__ == "__main__":
    main()

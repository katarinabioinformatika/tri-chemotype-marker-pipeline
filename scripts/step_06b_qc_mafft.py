#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys
from Bio import SeqIO

BASE_DIR = Path(__file__).resolve().parents[1]

# Use the original MAFFT headers so its "_R_" orientation markers are retained.
ALIGNMENT_FILE = (
    BASE_DIR
    / "results"
    / "mauve"
    / "core_TRI_cluster"
    / "tri_core_all_mafft.fasta"
)

METADATA_FILE = (
    BASE_DIR
    / "results"
    / "extracted_loci"
    / "extracted_loci_metadata.csv"
)

OUT_DIR = (
    BASE_DIR
    / "results"
    / "alignment_qc"
    / "mafft"
    / "core_TRI_cluster"
)

ALIGNMENT_METADATA_FILE = (
    OUT_DIR
    / "core_TRI_mafft_alignment_metadata.csv"
)

# Detect an incomplete or altered version of the current alignment.
EXPECTED_SEQUENCE_COUNT = 27

# Valid nucleotide, unresolved-base, and alignment-gap characters.
ALLOWED_CHARACTERS = set(
    "ACGTN-"
)

DISCOVERY_CHEMOTYPES = {
    "3-ADON",
    "15-ADON",
    "NIV",
    "T-2/HT-2",
}

def normalise_chemotype(value):
    """Standardise chemotype text imported from the metadata table."""

    if pd.isna(value):
        return ""

    return (
        str(value)
        .replace("\xa0", " ")
        .strip()
    )

def assign_analysis_role(chemotype):
    """Assign a downstream analysis role according to chemotype evidence."""

    if chemotype in DISCOVERY_CHEMOTYPES:
        return "discovery"

    if chemotype == "T-2/HT-2 *":
        return "validation_unconfirmed"

    if chemotype == "NP":
        return "negative_control"

    if chemotype == "Not reported":
        return "unlabelled"

    return "other"

def extract_genome_id(fasta_header):
    """Extract the genome accession without altering the aligned sequence."""

    genome_id = str(fasta_header).split("|")[0].strip()

    if genome_id.startswith("_R_"):
        genome_id = genome_id[3:]

    return genome_id

def was_reverse_complemented(fasta_header):
    """Determine whether MAFFT marked a sequence as reverse-complemented."""

    first_field = (
        str(fasta_header)
        .split("|")[0]
        .strip()
    )

    return first_field.startswith("_R_")

def main():
    """Check structural consistency of the core TRI-cluster MAFFT alignment.

    This does not assess whether MAFFT produced the optimal biological
    alignment.
    """

    if not ALIGNMENT_FILE.exists():
        sys.exit(
            f"ERROR: MAFFT alignment FASTA not found:\n{ALIGNMENT_FILE}"
        )

    if not METADATA_FILE.exists():
        sys.exit(
            f"ERROR: Extraction metadata file not found:\n{METADATA_FILE}"
        )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    records = list(
        SeqIO.parse(
            ALIGNMENT_FILE,
            "fasta"
        )
    )

    if not records:
        sys.exit(
            "ERROR: No FASTA records were found in the MAFFT alignment."
        )

    print("\nCore TRI-cluster MAFFT alignment sanity check")
    print("---------------------------------------------")

    sequence_count = len(records)

    print(
        f"FASTA records:                       {sequence_count}"
    )

    if sequence_count == EXPECTED_SEQUENCE_COUNT:
        print(
            f"Expected sequence count:             "
            f"PASS ({EXPECTED_SEQUENCE_COUNT})"
        )
    else:
        print(
            f"Expected sequence count:             "
            f"FAIL (expected {EXPECTED_SEQUENCE_COUNT})"
        )

    headers = [
        record.id
        for record in records
    ]

    unique_header_count = len(
        set(headers)
    )

    print(
        f"Unique FASTA headers:                "
        f"{unique_header_count}/{sequence_count}"
    )

    if unique_header_count != sequence_count:
        print(
            "WARNING: Duplicate FASTA headers were detected."
        )

    genome_ids = [
        extract_genome_id(
            record.id
        )
        for record in records
    ]

    unique_genome_count = len(
        set(genome_ids)
    )

    print(
        f"Unique Genome_ID values:             "
        f"{unique_genome_count}/{sequence_count}"
    )

    if unique_genome_count != sequence_count:
        print(
            "WARNING: At least one Genome_ID occurs more than once."
        )

    reverse_flags = [
        was_reverse_complemented(
            record.id
        )
        for record in records
    ]

    reverse_count = sum(
        reverse_flags
    )

    print(
        f"Reverse-complemented by MAFFT:        {reverse_count}"
    )

    if reverse_count > 0:

        print(
            "\nGenomes reverse-complemented by MAFFT:"
        )

        for record in records:

            if was_reverse_complemented(
                record.id
            ):
                print(
                    f"  {extract_genome_id(record.id)}"
                )

    alignment_lengths = [
        len(record.seq)
        for record in records
    ]

    unique_lengths = sorted(
        set(alignment_lengths)
    )

    if len(unique_lengths) == 1:

        alignment_length = unique_lengths[0]

        print(
            f"\nAlignment length:                    "
            f"{alignment_length:,} columns"
        )

        print(
            "All sequences same length:           PASS"
        )

    else:

        print(
            "\nAll sequences same length:           FAIL"
        )

        print(
            f"Observed lengths:                    {unique_lengths}"
        )

    invalid_character_records = []

    for record in records:

        sequence = str(
            record.seq
        ).upper()

        observed_characters = set(
            sequence
        )

        invalid_characters = (
            observed_characters
            - ALLOWED_CHARACTERS
        )

        if invalid_characters:

            invalid_character_records.append(
                (
                    record.id,
                    invalid_characters
                )
            )

    if not invalid_character_records:

        print(
            "Sequence characters:                 PASS"
        )

    else:

        print(
            "Sequence characters:                 FAIL"
        )

        for header, characters in invalid_character_records:

            print(
                f"  {header}: {sorted(characters)}"
            )

    metadata = pd.read_csv(
        METADATA_FILE
    )

    required_columns = [
        "Genome_ID",
        "Species",
        "Chemotype",
        "Chromosome_or_Contig",
        "Locus_class",
        "Genes",
        "Locus_start",
        "Locus_end",
        "Extract_start",
        "Extract_end",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in metadata.columns
    ]

    if missing_columns:

        sys.exit(
            "ERROR: Metadata table is missing required columns: "
            f"{missing_columns}"
        )

    core_metadata = metadata[
        metadata["Locus_class"] == "core_TRI_cluster"
    ].copy()

    print(
        f"Core TRI metadata records:           {len(core_metadata)}"
    )

    core_metadata["Genome_ID"] = (
        core_metadata["Genome_ID"]
        .astype(str)
        .str.strip()
    )

    core_metadata["Chemotype"] = (
        core_metadata["Chemotype"]
        .apply(normalise_chemotype)
    )

    duplicate_metadata_ids = (
        core_metadata[
            core_metadata["Genome_ID"].duplicated(
                keep=False
            )
        ]
    )

    if duplicate_metadata_ids.empty:

        print(
            "Core metadata Genome_IDs:            PASS (unique)"
        )

    else:

        print(
            "Core metadata Genome_IDs:            FAIL (duplicates found)"
        )

        print(
            duplicate_metadata_ids[
                [
                    "Genome_ID",
                    "Chromosome_or_Contig",
                    "Locus_start",
                    "Locus_end",
                ]
            ].to_string(
                index=False
            )
        )

    alignment_genome_set = set(
        genome_ids
    )

    metadata_genome_set = set(
        core_metadata["Genome_ID"]
    )

    missing_from_metadata = sorted(
        alignment_genome_set
        - metadata_genome_set
    )

    missing_from_alignment = sorted(
        metadata_genome_set
        - alignment_genome_set
    )

    matched_genomes = (
        alignment_genome_set
        & metadata_genome_set
    )

    print(
        f"Alignment genomes matched:            "
        f"{len(matched_genomes)}/{sequence_count}"
    )

    if missing_from_metadata:

        print(
            "\nGenome IDs present in MAFFT alignment "
            "but missing from metadata:"
        )

        for genome_id in missing_from_metadata:
            print(
                f"  {genome_id}"
            )

    else:

        print(
            "Missing from metadata:               0"
        )

    if missing_from_alignment:

        print(
            "\nCore TRI genomes present in metadata "
            "but missing from MAFFT alignment:"
        )

        for genome_id in missing_from_alignment:
            print(
                f"  {genome_id}"
            )

    else:

        print(
            "Missing from alignment:              0"
        )

    alignment_rows = []

    for position, record in enumerate(
        records,
        start=1
    ):

        genome_id = extract_genome_id(
            record.id
        )

        reverse_complemented = was_reverse_complemented(
            record.id
        )

        sequence = str(
            record.seq
        ).upper()

        gap_count = sequence.count("-")
        n_count = sequence.count("N")
        non_gap_length = len(sequence) - gap_count

        alignment_rows.append(
            {
                "Alignment_order":
                    position,

                "Genome_ID":
                    genome_id,

                "MAFFT_header":
                    record.id,

                "MAFFT_reverse_complemented":
                    reverse_complemented,

                "Alignment_length":
                    len(sequence),

                "Gap_count":
                    gap_count,

                "N_count":
                    n_count,

                "Non_gap_length":
                    non_gap_length,
            }
        )

    alignment_table = pd.DataFrame(
        alignment_rows
    )

    columns_to_attach = [
        "Genome_ID",
        "Species",
        "Chemotype",
        "Chromosome_or_Contig",
        "Genes",
        "Locus_start",
        "Locus_end",
        "Extract_start",
        "Extract_end",
    ]

    alignment_table = alignment_table.merge(
        core_metadata[
            columns_to_attach
        ],
        on="Genome_ID",
        how="left",
        validate="one_to_one"
    )

    alignment_table["Analysis_role"] = (
        alignment_table["Chemotype"]
        .apply(assign_analysis_role)
    )

    alignment_table["Gap_fraction"] = (
        alignment_table["Gap_count"]
        / alignment_table["Alignment_length"]
    )

    print(
        "\nChemotype distribution"
    )

    print(
        "-----------------------"
    )

    print(
        alignment_table[
            "Chemotype"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )

    print(
        "\nAnalysis roles"
    )

    print(
        "--------------"
    )

    print(
        alignment_table[
            "Analysis_role"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )

    print(
        "\nMAFFT alignment sequence mapping"
    )

    print(
        "--------------------------------"
    )

    display_columns = [
        "Alignment_order",
        "Genome_ID",
        "Species",
        "Chemotype",
        "Analysis_role",
        "MAFFT_reverse_complemented",
        "Non_gap_length",
        "Gap_fraction",
    ]

    print(
        alignment_table[
            display_columns
        ]
        .to_string(
            index=False
        )
    )

    alignment_table.to_csv(
        ALIGNMENT_METADATA_FILE,
        index=False
    )

    problems = []

    if sequence_count != EXPECTED_SEQUENCE_COUNT:

        problems.append(
            "unexpected number of FASTA records"
        )

    if unique_header_count != sequence_count:

        problems.append(
            "duplicate FASTA headers"
        )

    if unique_genome_count != sequence_count:

        problems.append(
            "duplicate Genome_ID values"
        )

    if len(unique_lengths) != 1:

        problems.append(
            "aligned sequences have different lengths"
        )

    if invalid_character_records:

        problems.append(
            "unexpected sequence characters"
        )

    if missing_from_metadata:

        problems.append(
            "alignment genomes missing from metadata"
        )

    if missing_from_alignment:

        problems.append(
            "core TRI metadata genomes missing from alignment"
        )

    if not duplicate_metadata_ids.empty:

        problems.append(
            "duplicate core TRI Genome_ID records in metadata"
        )

    print(
        "\n---------------------------------------------"
    )

    if problems:

        print(
            "SANITY CHECK RESULT: FAIL"
        )

        print(
            "\nProblems detected:"
        )

        for problem in problems:

            print(
                f"  - {problem}"
            )

    else:

        print(
            "SANITY CHECK RESULT: PASS"
        )

        print(
            "\nThe MAFFT alignment is structurally consistent "
            "and suitable for downstream comparison."
        )

    print(
        f"\nAlignment metadata saved to:\n"
        f"{ALIGNMENT_METADATA_FILE}"
    )

if __name__ == "__main__":
    main()

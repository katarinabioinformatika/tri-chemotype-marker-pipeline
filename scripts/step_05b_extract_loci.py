#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys
from Bio import SeqIO

BASE_DIR = Path(__file__).resolve().parents[1]

LOCI_FILE = (
    BASE_DIR
    / "results"
    / "classified_loci"
    / "tri_loci_classified.csv"
)

GENOME_DIR = BASE_DIR / "genomes"

OUT_DIR = BASE_DIR / "results" / "extracted_loci"

METADATA_FILE = OUT_DIR / "extracted_loci_metadata.csv"

# Include surrounding sequence for downstream structural and marker analysis.
FLANK_BP = 5_000

# Exclude heterogeneous "other_TRI_locus" regions from structural comparison.
CLASSES_TO_EXTRACT = {
    "core_TRI_cluster",
    "TRI101_locus",
    "TRI1_TRI16_locus",
}

def find_genome_file(genome_id):
    """Locate the genome FASTA file corresponding to one genome accession."""

    candidates = list(
        GENOME_DIR.glob(
            f"{genome_id}*.fna"
        )
    )

    if not candidates:
        candidates = list(
            GENOME_DIR.glob(
                f"{genome_id}*.fa"
            )
        )

    if not candidates:
        candidates = list(
            GENOME_DIR.glob(
                f"{genome_id}*.fasta"
            )
        )

    if not candidates:
        return None

    return candidates[0]

def load_genome_sequences(genome_file):
    """Load all chromosome or contig sequences from one genome FASTA file."""

    return SeqIO.to_dict(
        SeqIO.parse(
            genome_file,
            "fasta"
        )
    )

def safe_filename(text):
    """Convert text into a form that is safe to use inside output filenames."""

    return (
        str(text)
        .replace("|", "_")
        .replace("/", "_")
        .replace(":", "_")
        .replace(" ", "_")
    )

def main():
    """Extract classified TRI loci together with flanking genomic sequence."""

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    if not LOCI_FILE.exists():
        sys.exit(
            f"ERROR: Classified loci file not found: {LOCI_FILE}"
        )

    loci = pd.read_csv(
        LOCI_FILE
    )

    if loci.empty:
        sys.exit(
            f"ERROR: Classified loci file contains no records: {LOCI_FILE}"
        )

    required_cols = [
        "Genome_ID",
        "Chromosome_or_Contig",
        "Locus_ID",
        "Locus_start",
        "Locus_end",
        "Genes",
        "Locus_class",
    ]

    missing = [
        column
        for column in required_cols
        if column not in loci.columns
    ]

    if missing:
        sys.exit(
            f"ERROR: Missing required columns: {missing}"
        )

    loci = loci[
        loci["Locus_class"].isin(
            CLASSES_TO_EXTRACT
        )
    ].copy()

    if loci.empty:
        sys.exit(
            "ERROR: No loci matched the selected locus classes."
        )

    extracted_records = []

    genome_cache = {}

    for _, row in loci.iterrows():

        genome_id = str(
            row["Genome_ID"]
        )

        contig = str(
            row["Chromosome_or_Contig"]
        )

        locus_class = str(
            row["Locus_class"]
        )

        genome_file = find_genome_file(
            genome_id
        )

        if genome_file is None:
            print(
                "WARNING: Genome FASTA not found for "
                f"{genome_id}; skipping."
            )
            continue

        if genome_id not in genome_cache:
            genome_cache[genome_id] = (
                load_genome_sequences(
                    genome_file
                )
            )

        genome_seqs = genome_cache[
            genome_id
        ]

        if contig not in genome_seqs:
            print(
                f"WARNING: Contig {contig} not found in "
                f"{genome_file.name}; skipping."
            )
            continue

        seq_record = genome_seqs[
            contig
        ]

        contig_length = len(
            seq_record.seq
        )

        locus_start = int(
            row["Locus_start"]
        )

        locus_end = int(
            row["Locus_end"]
        )

        extract_start = max(
            1,
            locus_start - FLANK_BP
        )

        extract_end = min(
            contig_length,
            locus_end + FLANK_BP
        )

        extracted_seq = (
            seq_record.seq[
                extract_start - 1:
                extract_end
            ]
        )

        class_dir = (
            OUT_DIR
            / locus_class
        )

        class_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        output_name = (
            f"{safe_filename(genome_id)}_"
            f"{safe_filename(contig)}_"
            f"{safe_filename(locus_class)}_"
            f"{extract_start}_{extract_end}.fna"
        )

        output_path = (
            class_dir
            / output_name
        )

        header = (
            f"{genome_id}|"
            f"{contig}|"
            f"{locus_class}|"
            f"{extract_start}-{extract_end}|"
            f"locus={locus_start}-{locus_end}|"
            f"genes={row['Genes']}"
        )

        with open(
            output_path,
            "w"
        ) as out:

            out.write(
                f">{header}\n"
            )

            for i in range(
                0,
                len(extracted_seq),
                80
            ):
                out.write(
                    str(
                        extracted_seq[
                            i:i + 80
                        ]
                    )
                    + "\n"
                )

        extracted_records.append(
            {
                "Genome_ID":
                    genome_id,

                "Species":
                    row.get(
                        "Species",
                        ""
                    ),

                "Chemotype":
                    row.get(
                        "Chemotype",
                        ""
                    ),

                "Toxin_System":
                    row.get(
                        "Toxin_System",
                        ""
                    ),

                "Strain":
                    row.get(
                        "Strain",
                        ""
                    ),

                "Chromosome_or_Contig":
                    contig,

                "Locus_ID":
                    row["Locus_ID"],

                "Locus_class":
                    locus_class,

                "Genes":
                    row["Genes"],

                "Locus_start":
                    locus_start,

                "Locus_end":
                    locus_end,

                "Locus_length_bp":
                    row.get(
                        "Locus_length_bp",
                        ""
                    ),

                "Extract_start":
                    extract_start,

                "Extract_end":
                    extract_end,

                "Extract_length_bp":
                    (
                        extract_end
                        - extract_start
                        + 1
                    ),

                "Flank_bp":
                    FLANK_BP,

                "Genome_fasta":
                    str(
                        genome_file
                    ),

                "Output_fasta":
                    str(
                        output_path
                    ),
            }
        )

    metadata = pd.DataFrame(
        extracted_records
    )

    metadata.to_csv(
        METADATA_FILE,
        index=False
    )

    print("\nDone.")

    print(
        "Loci selected for extraction: "
        f"{len(loci)}"
    )

    print(
        "Loci successfully extracted: "
        f"{len(metadata)}"
    )

    print(
        f"Output directory: {OUT_DIR}"
    )

    print(
        f"Metadata saved to: {METADATA_FILE}"
    )

    if len(metadata) > 0:

        print(
            "\nExtracted locus classes:"
        )

        print(
            metadata[
                "Locus_class"
            ]
            .value_counts()
            .to_string()
        )

if __name__ == "__main__":
    main()

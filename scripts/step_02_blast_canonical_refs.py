#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys


# Project paths
BASE_DIR = Path(__file__).resolve().parents[1]
GENOME_DIR = BASE_DIR / "genomes"
REF_FASTA = BASE_DIR / "reference_genes" / "canonical_TRI_references.fna"
DB_DIR = BASE_DIR / "databases" / "genome_blastdb"
OUT_DIR = BASE_DIR / "results" / "tri_blast"


# BLAST filtering parameters
EVALUE = "1e-20"
PERC_IDENTITY = "50"


def run(cmd):
    """Run an external command and stop if it fails."""
    print("Running:", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def main():
    """Create genome BLAST databases and search for canonical gene references."""
    DB_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not REF_FASTA.exists():
        sys.exit(f"ERROR: reference FASTA not found: {REF_FASTA}")

    genomes = sorted(GENOME_DIR.glob("*.fna"))

    if not genomes:
        sys.exit(f"ERROR: no .fna genome files found in: {GENOME_DIR}")

    print(f"Found {len(genomes)} genomes")
    print(f"Reference FASTA: {REF_FASTA}")

    for genome in genomes:
        base = genome.stem
        db_prefix = DB_DIR / base

        print(f"\n=== Creating BLAST database for {base} ===")

        run([
            "makeblastdb",
            "-in", genome,
            "-dbtype", "nucl",
            "-out", db_prefix,
        ])

        out_file = OUT_DIR / f"{base}_tri_hits.tsv"

        print(f"=== BLAST canonical references against {base} ===")

        run([
            "blastn",
            "-query", REF_FASTA,
            "-db", db_prefix,
            "-out", out_file,
            "-evalue", EVALUE,
            "-perc_identity", PERC_IDENTITY,
            "-outfmt",
            (
                "6 qseqid sseqid pident length qlen "
                "qstart qend sstart send evalue bitscore"
            ),
        ])

    print("\nDone.")
    print(f"Results saved in: {OUT_DIR}")


if __name__ == "__main__":
    main()
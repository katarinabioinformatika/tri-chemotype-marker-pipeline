#!/usr/bin/env python3

from pathlib import Path
import gzip
import shutil

import pandas as pd
import requests


# Project paths
PROJECT = Path(__file__).resolve().parents[1]
METADATA = PROJECT / "metadata" / "metadata.csv"
GENOMES = PROJECT / "genomes"

GENOMES.mkdir(exist_ok=True)


# NCBI assembly summaries
SUMMARY_URLS = [
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt",
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt",
]


def load_summaries():
    """Return assembly-accession-to-FTP-path mappings from NCBI summaries."""
    records = {}

    for url in SUMMARY_URLS:
        print(f"Reading {url}")
        response = requests.get(url, timeout=120)
        response.raise_for_status()

        for line in response.text.splitlines():
            if line.startswith("#"):
                continue

            cols = line.split("\t")
            accession = cols[0].strip()
            ftp_path = cols[19].strip().rstrip("/")

            if ftp_path and ftp_path != "na":
                records[accession] = ftp_path

    return records


def download_and_unzip(url, gz_path, output_path):
    """Download and decompress a gzip file, returning whether it succeeded."""
    if output_path.exists():
        print(f"Already exists: {output_path.name}")
        return True

    print(f"Downloading: {url}")
    response = requests.get(url, stream=True, timeout=300)

    if response.status_code == 404:
        print(f"Not available: {url}")
        return False

    response.raise_for_status()

    with open(gz_path, "wb") as out:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                out.write(chunk)

    with gzip.open(gz_path, "rb") as f_in:
        with open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()

    print(f"Saved: {output_path.name}")
    return True


def download_genome_file(accession, ftp_path):
    """Download the genomic FASTA file for one NCBI assembly."""
    ftp_path = ftp_path.rstrip("/")
    folder_name = ftp_path.split("/")[-1]

    fasta_url = f"{ftp_path}/{folder_name}_genomic.fna.gz"
    fasta_gz = GENOMES / f"{accession}_genomic.fna.gz"
    fasta_out = GENOMES / f"{accession}_genomic.fna"

    print(f"\nProcessing {accession}")

    success = download_and_unzip(fasta_url, fasta_gz, fasta_out)

    if not success:
        print(f"Genome download failed for {accession}")


def main():
    """Download all genome assemblies listed in the metadata table."""
    df = pd.read_csv(METADATA)

    if "Genome_ID" not in df.columns:
        raise ValueError(
            "metadata.csv must contain a column called Genome_ID"
        )

    records = load_summaries()

    accessions = (
        df["Genome_ID"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    for accession in accessions:
        if not accession:
            continue

        if accession not in records:
            print(f"Not found in NCBI assembly summaries: {accession}")
            continue

        download_genome_file(accession, records[accession])

    print("\nDone.")


if __name__ == "__main__":
    main()
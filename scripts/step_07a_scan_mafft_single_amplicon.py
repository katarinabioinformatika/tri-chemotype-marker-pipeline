#!/usr/bin/env python3

"""Scan a MAFFT TRI-core alignment for chemotype-separating amplicons."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE_DIR = Path(__file__).resolve().parents[1]

DEFAULT_ALIGNMENT = (
    BASE_DIR
    / "results"
    / "mauve"
    / "core_TRI_cluster"
    / "tri_core_all_mafft.fasta"
)

DEFAULT_METADATA = (
    BASE_DIR
    / "results"
    / "extracted_loci"
    / "extracted_loci_metadata.csv"
)

DEFAULT_OUTPUT = (
    BASE_DIR
    / "results"
    / "chemotype_separation"
    / "diagnostic_amplicon_scan"
)

# Target physical lengths count A/C/G/T bases rather than alignment columns.
TARGET_LENGTHS_BP = (
    1500,
    2000,
    2500,
    3000,
    3500,
    4000,
    4500,
    5000,
)

MIN_INDIVIDUAL_LENGTH_FRACTION = 0.80

# Candidate starts are tested every 100 alignment columns.
SCAN_STEP_COLUMNS = 100

# Broad terminal regions used to assess primer-site feasibility.
END_REGION_COLUMNS = 100

MIN_CANDIDATE_COVERAGE = 0.85
MIN_END_COVERAGE = 0.90

MIN_END_CONSERVATION = 0.90

# Require this proportion of the shorter sequence to be comparable.
MIN_PAIR_COMPARABLE_FRACTION = 0.70

# Diagnostic sites require sufficient coverage and within-group consensus.
MIN_DIAGNOSTIC_GROUP_COVERAGE = 0.80
MIN_DIAGNOSTIC_CONSENSUS = 0.75

TOP_CANDIDATES_TO_EXPORT = 10

GENOME_ID_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)", re.IGNORECASE)
VALID_BASES = set("ACGT")

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan a MAFFT TRI-core alignment for chemotype-separating amplicons."
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        default=DEFAULT_ALIGNMENT,
        help=f"Aligned FASTA file (default: {DEFAULT_ALIGNMENT})",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=DEFAULT_METADATA,
        help=f"Metadata CSV containing Genome_ID, Species and Chemotype (default: {DEFAULT_METADATA})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output directory (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=TOP_CANDIDATES_TO_EXPORT,
        help=f"Number of top candidates to export (default: {TOP_CANDIDATES_TO_EXPORT})",
    )
    return parser.parse_args()

def extract_genome_id(text: str) -> str | None:
    """Recover a GCA/GCF accession from a FASTA header or other text."""
    match = GENOME_ID_PATTERN.search(str(text))
    return match.group(1).upper() if match else None

def normalise_chemotype(value: object) -> str | None:
    """Convert metadata labels into four computational chemotype groups.

    Confirmed and literature-only T-2/HT-2 labels map to the same group;
    their evidence status remains available separately in the output.
    """
    if pd.isna(value):
        return None

    text = str(value).strip().upper()
    compact = re.sub(r"[\s_*]", "", text)
    compact = compact.replace("–", "-").replace("—", "-")

    if compact in {"3-ADON", "3ADON"}:
        return "3-ADON"
    if compact in {"15-ADON", "15ADON"}:
        return "15-ADON"
    if compact == "NIV":
        return "NIV"
    if "T-2/HT-2" in compact or "T2/HT2" in compact:
        return "T-2/HT-2"
    return None

def safe_name(text: object) -> str:
    """Create a compact text value suitable for FASTA headers and filenames."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text).strip()).strip("_")

def load_alignment(path: Path) -> list[SeqRecord]:
    if not path.exists():
        sys.exit(f"ERROR: MAFFT alignment not found: {path}")

    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        sys.exit(f"ERROR: No FASTA records were found in: {path}")

    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        sys.exit(
            "ERROR: The FASTA records do not all have the same length. "
            "The input must be a multiple-sequence alignment."
        )

    genome_ids = [extract_genome_id(record.id) for record in records]
    if any(genome_id is None for genome_id in genome_ids):
        bad = [record.id for record, genome_id in zip(records, genome_ids) if genome_id is None]
        sys.exit(f"ERROR: Could not recover Genome_ID from FASTA headers: {bad}")

    duplicates = sorted({x for x in genome_ids if genome_ids.count(x) > 1})
    if duplicates:
        sys.exit(f"ERROR: Duplicate Genome_ID values in the alignment: {duplicates}")

    return records

def load_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: Metadata CSV not found: {path}")

    metadata = pd.read_csv(path, encoding="utf-8-sig")
    required = {"Genome_ID", "Species", "Chemotype"}
    missing = sorted(required.difference(metadata.columns))
    if missing:
        sys.exit(f"ERROR: Metadata is missing required columns: {missing}")

    metadata = metadata.copy()
    metadata["Genome_ID"] = metadata["Genome_ID"].astype(str).str.strip().str.upper()

    if "Locus_class" in metadata.columns:
        core = metadata[metadata["Locus_class"].eq("core_TRI_cluster")].copy()
        if not core.empty:
            metadata = core

    duplicate_ids = metadata.loc[metadata["Genome_ID"].duplicated(), "Genome_ID"].unique()
    if len(duplicate_ids):
        conflicts = []
        for genome_id in duplicate_ids:
            subset = metadata.loc[metadata["Genome_ID"].eq(genome_id), ["Species", "Chemotype"]]
            if len(subset.drop_duplicates()) > 1:
                conflicts.append(genome_id)
        if conflicts:
            sys.exit(f"ERROR: Conflicting metadata rows for Genome_ID values: {conflicts}")
        metadata = metadata.drop_duplicates(subset="Genome_ID", keep="first")

    metadata["Separation_group"] = metadata["Chemotype"].map(normalise_chemotype)
    metadata["Evidence_status"] = np.where(
        metadata["Chemotype"].astype(str).str.contains(r"\*", regex=True),
        "inferred",
        "confirmed",
    )
    return metadata

def build_discovery_table(
    records: list[SeqRecord], metadata: pd.DataFrame
) -> tuple[list[SeqRecord], pd.DataFrame]:
    """Join aligned sequences to metadata and retain the four discovery groups."""
    rows = []
    record_by_id = {}

    metadata_by_id = metadata.set_index("Genome_ID", drop=False)
    for record in records:
        genome_id = extract_genome_id(record.id)
        record_by_id[genome_id] = record
        if genome_id not in metadata_by_id.index:
            sys.exit(f"ERROR: No metadata match for aligned genome: {genome_id}")

        row = metadata_by_id.loc[genome_id].to_dict()
        row["FASTA_header"] = record.description
        rows.append(row)

    joined = pd.DataFrame(rows)
    discovery = joined[joined["Separation_group"].notna()].copy()
    discovery_records = [record_by_id[x] for x in discovery["Genome_ID"]]
    discovery = discovery.reset_index(drop=True)

    required_groups = {"3-ADON", "15-ADON", "NIV", "T-2/HT-2"}
    observed_groups = set(discovery["Separation_group"])
    missing_groups = sorted(required_groups.difference(observed_groups))
    if missing_groups:
        sys.exit(f"ERROR: No discovery genomes were found for groups: {missing_groups}")

    return discovery_records, discovery

def encode_alignment(records: list[SeqRecord]) -> np.ndarray:
    """Encode A/C/G/T as 0/1/2/3 and gaps or ambiguous symbols as 4."""
    lookup = np.full(256, 4, dtype=np.uint8)
    for code, base in enumerate(b"ACGT"):
        lookup[base] = code
        lookup[ord(chr(base).lower())] = code

    matrix = np.empty((len(records), len(records[0].seq)), dtype=np.uint8)
    for row, record in enumerate(records):
        raw = np.frombuffer(str(record.seq).encode("ascii"), dtype=np.uint8)
        matrix[row] = lookup[raw]
    return matrix

def cumulative_counts(boolean_matrix: np.ndarray) -> np.ndarray:
    """Return row-wise prefix sums with an initial zero column."""
    prefix = np.zeros(
        (boolean_matrix.shape[0], boolean_matrix.shape[1] + 1),
        dtype=np.int32,
    )
    np.cumsum(boolean_matrix, axis=1, dtype=np.int32, out=prefix[:, 1:])
    return prefix

def make_pair_prefixes(encoded: np.ndarray):
    """Precalculate comparable-base and mismatch counts for all genome pairs."""
    pair_i, pair_j = np.triu_indices(encoded.shape[0], k=1)
    comparable = (encoded[pair_i] < 4) & (encoded[pair_j] < 4)
    mismatch = comparable & (encoded[pair_i] != encoded[pair_j])
    return pair_i, pair_j, cumulative_counts(comparable), cumulative_counts(mismatch)

def column_statistics(encoded: np.ndarray):
    """Calculate base coverage and consensus conservation for every column."""
    n_genomes = encoded.shape[0]
    counts = np.stack([(encoded == base).sum(axis=0) for base in range(4)])
    valid = counts.sum(axis=0)
    coverage = valid / n_genomes
    consensus = np.divide(
        counts.max(axis=0),
        valid,
        out=np.zeros(valid.shape, dtype=float),
        where=valid > 0,
    )

    coverage_prefix = np.concatenate(([0.0], np.cumsum(coverage)))
    consensus_prefix = np.concatenate(([0.0], np.cumsum(consensus)))
    return coverage_prefix, consensus_prefix

def interval_mean(prefix: np.ndarray, start: int, end: int) -> float:
    if end <= start:
        return float("nan")
    return float((prefix[end] - prefix[start]) / (end - start))

def median_physical_length(base_prefix: np.ndarray, start: int, end: int) -> float:
    return float(np.median(base_prefix[:, end] - base_prefix[:, start]))

def find_end_for_target(
    base_prefix: np.ndarray,
    start: int,
    target_bp: int,
    alignment_length: int,
) -> int | None:
    """Find the earliest end coordinate reaching a target median physical length."""
    if median_physical_length(base_prefix, start, alignment_length) < target_bp:
        return None

    low = start + 1
    high = alignment_length
    while low < high:
        middle = (low + high) // 2
        if median_physical_length(base_prefix, start, middle) >= target_bp:
            high = middle
        else:
            low = middle + 1
    return low

def group_pair_mask(
    groups: np.ndarray,
    pair_i: np.ndarray,
    pair_j: np.ndarray,
    group_a: set[str],
    group_b: set[str] | None = None,
) -> np.ndarray:
    in_a_i = np.isin(groups[pair_i], list(group_a))
    in_a_j = np.isin(groups[pair_j], list(group_a))
    if group_b is None:
        return in_a_i & in_a_j

    in_b_i = np.isin(groups[pair_i], list(group_b))
    in_b_j = np.isin(groups[pair_j], list(group_b))
    return (in_a_i & in_b_j) | (in_b_i & in_a_j)

def calculate_margin(
    distances: np.ndarray,
    within_a: np.ndarray,
    within_b: np.ndarray,
    between: np.ndarray,
) -> tuple[float, float, float, float]:
    """Return max within A, max within B, min between, and their margin."""
    valid_a = distances[within_a & np.isfinite(distances)]
    valid_b = distances[within_b & np.isfinite(distances)]
    valid_between = distances[between & np.isfinite(distances)]
    if not len(valid_a) or not len(valid_b) or not len(valid_between):
        return (float("nan"),) * 4

    max_a = float(valid_a.max())
    max_b = float(valid_b.max())
    min_between = float(valid_between.min())
    margin = min_between - max(max_a, max_b)
    return max_a, max_b, min_between, margin

def find_diagnostic_sites(
    encoded: np.ndarray,
    groups: np.ndarray,
    start: int,
    end: int,
    group_a: set[str],
    group_b: set[str],
    group_a_name: str,
    group_b_name: str,
) -> pd.DataFrame:
    """
    Identify conserved, contrasting nucleotide positions between two groups.

    Gaps and ambiguous characters are excluded when calculating each group's
    consensus. Coverage is still checked against the full number of genomes in
    the group, preventing gap-rich positions from being called diagnostic.
    """
    mask_a = np.isin(groups, list(group_a))
    mask_b = np.isin(groups, list(group_b))
    sequences_a = encoded[mask_a, start:end]
    sequences_b = encoded[mask_b, start:end]

    counts_a = np.stack([(sequences_a == base).sum(axis=0) for base in range(4)])
    counts_b = np.stack([(sequences_b == base).sum(axis=0) for base in range(4)])
    valid_a = counts_a.sum(axis=0)
    valid_b = counts_b.sum(axis=0)

    coverage_a = valid_a / sequences_a.shape[0]
    coverage_b = valid_b / sequences_b.shape[0]
    consensus_code_a = counts_a.argmax(axis=0)
    consensus_code_b = counts_b.argmax(axis=0)
    consensus_fraction_a = np.divide(
        counts_a.max(axis=0),
        valid_a,
        out=np.zeros(valid_a.shape, dtype=float),
        where=valid_a > 0,
    )
    consensus_fraction_b = np.divide(
        counts_b.max(axis=0),
        valid_b,
        out=np.zeros(valid_b.shape, dtype=float),
        where=valid_b > 0,
    )

    diagnostic = (
        (coverage_a >= MIN_DIAGNOSTIC_GROUP_COVERAGE)
        & (coverage_b >= MIN_DIAGNOSTIC_GROUP_COVERAGE)
        & (consensus_fraction_a >= MIN_DIAGNOSTIC_CONSENSUS)
        & (consensus_fraction_b >= MIN_DIAGNOSTIC_CONSENSUS)
        & (consensus_code_a != consensus_code_b)
    )

    offsets = np.flatnonzero(diagnostic)
    bases = np.asarray(list("ACGT"))
    return pd.DataFrame(
        {
            "alignment_position": start + offsets + 1,
            "group_a": group_a_name,
            "group_b": group_b_name,
            "group_a_consensus": bases[consensus_code_a[offsets]],
            "group_b_consensus": bases[consensus_code_b[offsets]],
            "group_a_coverage": coverage_a[offsets],
            "group_b_coverage": coverage_b[offsets],
            "group_a_consensus_fraction": consensus_fraction_a[offsets],
            "group_b_consensus_fraction": consensus_fraction_b[offsets],
            "group_a_consensus_count": counts_a.max(axis=0)[offsets],
            "group_b_consensus_count": counts_b.max(axis=0)[offsets],
            "group_a_available_count": valid_a[offsets],
            "group_b_available_count": valid_b[offsets],
        }
    )

def candidate_diagnostic_sites(
    encoded: np.ndarray,
    groups: np.ndarray,
    start: int,
    end: int,
) -> dict[str, pd.DataFrame]:
    """Return diagnostic-position tables for the three planned contrasts."""
    return {
        "adon": find_diagnostic_sites(
            encoded,
            groups,
            start,
            end,
            {"3-ADON"},
            {"15-ADON"},
            "3-ADON",
            "15-ADON",
        ),
        "don_niv": find_diagnostic_sites(
            encoded,
            groups,
            start,
            end,
            {"3-ADON", "15-ADON"},
            {"NIV"},
            "DON",
            "NIV",
        ),
        "t2_typeb": find_diagnostic_sites(
            encoded,
            groups,
            start,
            end,
            {"T-2/HT-2"},
            {"3-ADON", "15-ADON", "NIV"},
            "T-2/HT-2",
            "TypeB",
        ),
    }

def scan_candidates(
    encoded: np.ndarray,
    discovery: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    n_genomes, alignment_length = encoded.shape
    valid = encoded < 4
    base_prefix = cumulative_counts(valid)
    pair_i, pair_j, comparable_prefix, mismatch_prefix = make_pair_prefixes(encoded)
    coverage_prefix, consensus_prefix = column_statistics(encoded)

    groups = discovery["Separation_group"].to_numpy()
    group_sets = {
        "3": {"3-ADON"},
        "15": {"15-ADON"},
        "DON": {"3-ADON", "15-ADON"},
        "NIV": {"NIV"},
        "T2": {"T-2/HT-2"},
        "TypeB": {"3-ADON", "15-ADON", "NIV"},
    }

    masks = {
        "within_3": group_pair_mask(groups, pair_i, pair_j, group_sets["3"]),
        "within_15": group_pair_mask(groups, pair_i, pair_j, group_sets["15"]),
        "between_adon": group_pair_mask(
            groups, pair_i, pair_j, group_sets["3"], group_sets["15"]
        ),
        "within_don": group_pair_mask(groups, pair_i, pair_j, group_sets["DON"]),
        "within_niv": group_pair_mask(groups, pair_i, pair_j, group_sets["NIV"]),
        "between_don_niv": group_pair_mask(
            groups, pair_i, pair_j, group_sets["DON"], group_sets["NIV"]
        ),
        "within_t2": group_pair_mask(groups, pair_i, pair_j, group_sets["T2"]),
        "within_typeb": group_pair_mask(groups, pair_i, pair_j, group_sets["TypeB"]),
        "between_t2_typeb": group_pair_mask(
            groups, pair_i, pair_j, group_sets["T2"], group_sets["TypeB"]
        ),
    }

    rows = []
    seen_intervals = set()

    for start in range(0, alignment_length, SCAN_STEP_COLUMNS):
        for target_bp in TARGET_LENGTHS_BP:
            end = find_end_for_target(base_prefix, start, target_bp, alignment_length)
            if end is None or (start, end) in seen_intervals:
                continue
            seen_intervals.add((start, end))

            physical_lengths = base_prefix[:, end] - base_prefix[:, start]
            median_length = float(np.median(physical_lengths))

            if (
                median_length < min(TARGET_LENGTHS_BP)
                or median_length > max(TARGET_LENGTHS_BP) + 25
            ):
                continue

            minimum_required_length = median_length * MIN_INDIVIDUAL_LENGTH_FRACTION

            if physical_lengths.min() < minimum_required_length:
                continue

            candidate_coverage = interval_mean(coverage_prefix, start, end)
            if candidate_coverage < MIN_CANDIDATE_COVERAGE:
                continue

            left_end = min(start + END_REGION_COLUMNS, end)
            right_start = max(start, end - END_REGION_COLUMNS)
            left_coverage = interval_mean(coverage_prefix, start, left_end)
            right_coverage = interval_mean(coverage_prefix, right_start, end)
            left_conservation = interval_mean(consensus_prefix, start, left_end)
            right_conservation = interval_mean(consensus_prefix, right_start, end)

            if min(left_coverage, right_coverage) < MIN_END_COVERAGE:
                continue
            if min(left_conservation, right_conservation) < MIN_END_CONSERVATION:
                continue

            comparable = comparable_prefix[:, end] - comparable_prefix[:, start]
            mismatches = mismatch_prefix[:, end] - mismatch_prefix[:, start]
            shorter_pair_length = np.minimum(physical_lengths[pair_i], physical_lengths[pair_j])
            minimum_comparable = np.maximum(
                1,
                np.ceil(shorter_pair_length * MIN_PAIR_COMPARABLE_FRACTION).astype(int),
            )
            distances = np.divide(
                mismatches,
                comparable,
                out=np.full(comparable.shape, np.nan, dtype=float),
                where=comparable >= minimum_comparable,
            )

            adon = calculate_margin(
                distances,
                masks["within_3"],
                masks["within_15"],
                masks["between_adon"],
            )
            don_niv = calculate_margin(
                distances,
                masks["within_don"],
                masks["within_niv"],
                masks["between_don_niv"],
            )
            t2_typeb = calculate_margin(
                distances,
                masks["within_t2"],
                masks["within_typeb"],
                masks["between_t2_typeb"],
            )

            margins = [adon[3], don_niv[3], t2_typeb[3]]
            balanced = float(np.min(margins)) if np.all(np.isfinite(margins)) else np.nan

            diagnostic_tables = candidate_diagnostic_sites(
                encoded,
                groups,
                start,
                end,
            )
            diagnostic_adon = len(diagnostic_tables["adon"])
            diagnostic_don_niv = len(diagnostic_tables["don_niv"])
            diagnostic_t2_typeb = len(diagnostic_tables["t2_typeb"])
            diagnostic_counts = [
                diagnostic_adon,
                diagnostic_don_niv,
                diagnostic_t2_typeb,
            ]
            balanced_diagnostic_sites = min(diagnostic_counts)
            total_diagnostic_sites = sum(diagnostic_counts)
            median_length_kb = median_length / 1000.0
            diagnostic_density_adon = diagnostic_adon / median_length_kb
            diagnostic_density_don_niv = diagnostic_don_niv / median_length_kb
            diagnostic_density_t2_typeb = diagnostic_t2_typeb / median_length_kb
            balanced_diagnostic_density = min(
                diagnostic_density_adon,
                diagnostic_density_don_niv,
                diagnostic_density_t2_typeb,
            )
            positive_margin_count = int(
                np.sum(np.asarray(margins, dtype=float) > 0)
            )

            rows.append(
                {
                    "alignment_start": start + 1,
                    "alignment_end": end,
                    "alignment_width_columns": end - start,
                    "target_length_bp": target_bp,
                    "median_length_bp": median_length,
                    "min_length_bp": int(physical_lengths.min()),
                    "max_length_bp": int(physical_lengths.max()),
                    "q10_length_bp": float(np.quantile(physical_lengths, 0.10)),
                    "q90_length_bp": float(np.quantile(physical_lengths, 0.90)),
                    "mean_candidate_coverage": candidate_coverage,
                    "left_end_coverage": left_coverage,
                    "right_end_coverage": right_coverage,
                    "left_end_conservation": left_conservation,
                    "right_end_conservation": right_conservation,
                    "diagnostic_sites_adon": diagnostic_adon,
                    "diagnostic_sites_don_niv": diagnostic_don_niv,
                    "diagnostic_sites_t2_typeb": diagnostic_t2_typeb,
                    "balanced_diagnostic_sites": balanced_diagnostic_sites,
                    "total_diagnostic_sites": total_diagnostic_sites,
                    "diagnostic_density_adon_per_kb": diagnostic_density_adon,
                    "diagnostic_density_don_niv_per_kb": diagnostic_density_don_niv,
                    "diagnostic_density_t2_typeb_per_kb": diagnostic_density_t2_typeb,
                    "balanced_diagnostic_density_per_kb": balanced_diagnostic_density,
                    "all_contrasts_have_diagnostic_sites": bool(
                        balanced_diagnostic_sites > 0
                    ),
                    "within_3adon_max": adon[0],
                    "within_15adon_max": adon[1],
                    "between_3adon_15adon_min": adon[2],
                    "margin_adon": adon[3],
                    "within_don_max": don_niv[0],
                    "within_niv_max": don_niv[1],
                    "between_don_niv_min": don_niv[2],
                    "margin_don_niv": don_niv[3],
                    "within_t2_max": t2_typeb[0],
                    "within_typeb_max": t2_typeb[1],
                    "between_t2_typeb_min": t2_typeb[2],
                    "margin_t2_typeb": t2_typeb[3],
                    "positive_margin_count": positive_margin_count,
                    "balanced_margin": balanced,
                    "all_margins_positive": bool(
                        np.all(np.asarray(margins, dtype=float) > 0)
                    ),
                    "minimum_pair_comparable_bp": int(comparable.min()),
                    "mean_pair_comparable_bp": float(comparable.mean()),
                }
            )

    candidates = pd.DataFrame(rows)
    context = {
        "base_prefix": base_prefix,
        "pair_i": pair_i,
        "pair_j": pair_j,
        "comparable_prefix": comparable_prefix,
        "mismatch_prefix": mismatch_prefix,
    }
    return candidates, context

def remove_highly_overlapping_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Retain the best representative when candidate intervals overlap by >=80%."""
    selected = []
    for index, row in candidates.iterrows():
        start = int(row["alignment_start"])
        end = int(row["alignment_end"])
        width = end - start + 1
        redundant = False
        for chosen in selected:
            overlap = max(0, min(end, chosen[1]) - max(start, chosen[0]) + 1)
            if overlap / min(width, chosen[2]) >= 0.80:
                redundant = True
                break
        if not redundant:
            selected.append((start, end, width, index))
    return candidates.loc[[x[3] for x in selected]].copy()

def distance_matrix_for_candidate(
    row: pd.Series,
    discovery: pd.DataFrame,
    context: dict,
) -> pd.DataFrame:
    start = int(row["alignment_start"]) - 1
    end = int(row["alignment_end"])
    pair_i = context["pair_i"]
    pair_j = context["pair_j"]
    comparable = context["comparable_prefix"][:, end] - context["comparable_prefix"][:, start]
    mismatches = context["mismatch_prefix"][:, end] - context["mismatch_prefix"][:, start]
    distances = np.divide(
        mismatches,
        comparable,
        out=np.full(comparable.shape, np.nan, dtype=float),
        where=comparable > 0,
    )

    labels = [
        f"{group}|{genome_id}"
        for group, genome_id in zip(
            discovery["Separation_group"], discovery["Genome_ID"]
        )
    ]
    matrix = np.zeros((len(labels), len(labels)), dtype=float)
    matrix[pair_i, pair_j] = distances
    matrix[pair_j, pair_i] = distances
    return pd.DataFrame(matrix, index=labels, columns=labels)

def export_candidate_fastas(
    row: pd.Series,
    rank: int,
    records: list[SeqRecord],
    discovery: pd.DataFrame,
    output_dir: Path,
) -> None:
    start = int(row["alignment_start"]) - 1
    end = int(row["alignment_end"])
    candidate_id = f"candidate_{rank:02d}_{start + 1}-{end}"
    candidate_dir = output_dir / "top_candidates" / candidate_id
    candidate_dir.mkdir(parents=True, exist_ok=True)

    aligned_records = []
    ungapped_records = []
    for record, (_, metadata_row) in zip(records, discovery.iterrows()):
        aligned_sequence = str(record.seq[start:end]).upper()
        ungapped_sequence = "".join(base for base in aligned_sequence if base in VALID_BASES)
        species = safe_name(metadata_row["Species"])
        group = safe_name(metadata_row["Separation_group"])
        evidence = safe_name(metadata_row["Evidence_status"])
        genome_id = metadata_row["Genome_ID"]
        header = f"{group}|{species}|{genome_id}|{evidence}"

        aligned_records.append(SeqRecord(Seq(aligned_sequence), id=header, description=""))
        ungapped_records.append(SeqRecord(Seq(ungapped_sequence), id=header, description=""))

    SeqIO.write(aligned_records, candidate_dir / f"{candidate_id}_aligned.fasta", "fasta")
    SeqIO.write(
        ungapped_records,
        candidate_dir / f"{candidate_id}_ungapped_for_Geneious.fasta",
        "fasta",
    )

def export_diagnostic_positions(
    row: pd.Series,
    encoded: np.ndarray,
    discovery: pd.DataFrame,
    candidate_dir: Path,
    candidate_id: str,
) -> None:
    """Export the actual diagnostic alignment positions for a top candidate."""
    start = int(row["alignment_start"]) - 1
    end = int(row["alignment_end"])
    groups = discovery["Separation_group"].to_numpy()
    tables = candidate_diagnostic_sites(encoded, groups, start, end)

    labelled_tables = []
    for contrast, table in tables.items():
        table = table.copy()
        table.insert(0, "contrast", contrast)
        labelled_tables.append(table)

    combined = pd.concat(labelled_tables, ignore_index=True)
    combined.to_csv(
        candidate_dir / f"{candidate_id}_diagnostic_positions.csv",
        index=False,
    )

def write_summary(
    output_dir: Path,
    alignment: Path,
    metadata_path: Path,
    discovery: pd.DataFrame,
    all_candidates: pd.DataFrame,
    shortlisted: pd.DataFrame,
) -> None:
    counts = discovery.groupby(["Separation_group", "Evidence_status"]).size()
    lines = [
        "Simplified chemotype-separating amplicon scan",
        "================================================",
        "",
        f"MAFFT alignment: {alignment}",
        f"Metadata:        {metadata_path}",
        f"Discovery genomes: {len(discovery)}",
        "",
        "Discovery groups:",
    ]
    for (group, evidence), count in counts.items():
        lines.append(f"  {group:12s} {evidence:9s}: {count}")

    diagnostic_all_three = (
        int(all_candidates["all_contrasts_have_diagnostic_sites"].sum())
        if len(all_candidates)
        else 0
    )
    positive = int(all_candidates["all_margins_positive"].sum()) if len(all_candidates) else 0
    lines.extend(
        [
            "",
            f"Candidates passing coverage/end filters: {len(all_candidates)}",
            f"Candidates with diagnostic sites in all three contrasts: {diagnostic_all_three}",
            f"Non-redundant shortlisted candidates:   {len(shortlisted)}",
            "",
            "Primary ranking:",
            "  Candidates are ranked first by the smallest diagnostic-site count",
            "  across ADON, DON/NIV and T-2/HT-2 versus Type B, followed by",
            "  diagnostic-site density and total diagnostic-site count.",
            "",
            "Supporting whole-sequence margin results:",
            f"  Candidates with all three margins > 0: {positive}",
            "  > 0  complete separation in the analysed discovery genomes",
            "  = 0  groups touch",
            "  < 0  within-group and between-group distances overlap",
            "",
            "Distance margins are retained as supporting evidence and are not a",
            "pass/fail requirement. Exact primers must be designed and tested",
            "separately in Geneious Prime.",
        ]
    )
    (output_dir / "analysis_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

def main() -> None:
    args = parse_arguments()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_alignment(args.alignment.resolve())
    metadata = load_metadata(args.metadata.resolve())
    discovery_records, discovery = build_discovery_table(records, metadata)

    print("Chemotype-separating amplicon scan")
    print("-----------------------------------")
    print(f"Alignment records:   {len(records)}")
    print(f"Alignment length:    {len(records[0].seq):,} columns")
    print(f"Discovery genomes:   {len(discovery_records)}")
    print()
    print(discovery["Separation_group"].value_counts().sort_index().to_string())
    print()
    print("Scanning candidate regions...")

    discovery.to_csv(output_dir / "discovery_genomes_used.csv", index=False)
    encoded = encode_alignment(discovery_records)
    candidates, context = scan_candidates(encoded, discovery)

    if candidates.empty:
        write_summary(
            output_dir,
            args.alignment.resolve(),
            args.metadata.resolve(),
            discovery,
            candidates,
            candidates,
        )
        sys.exit(
            "No candidates passed the current coverage and conserved-end filters. "
            "Review analysis_summary.txt and relax the configurable thresholds if justified."
        )

    candidates = candidates.sort_values(
        by=[
            "balanced_diagnostic_sites",
            "balanced_diagnostic_density_per_kb",
            "total_diagnostic_sites",
            "positive_margin_count",
            "balanced_margin",
            "left_end_conservation",
            "right_end_conservation",
            "mean_candidate_coverage",
            "median_length_bp",
        ],
        ascending=[False, False, False, False, False, False, False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    candidates.insert(0, "overall_rank", np.arange(1, len(candidates) + 1))
    candidates.to_csv(output_dir / "all_ranked_candidates.csv", index=False)

    shortlisted = remove_highly_overlapping_candidates(candidates)
    shortlisted = shortlisted.head(max(1, args.top)).reset_index(drop=True)
    shortlisted.insert(0, "shortlist_rank", np.arange(1, len(shortlisted) + 1))
    shortlisted.to_csv(output_dir / "top_nonredundant_candidates.csv", index=False)

    for _, row in shortlisted.iterrows():
        rank = int(row["shortlist_rank"])
        export_candidate_fastas(
            row,
            rank,
            discovery_records,
            discovery,
            output_dir,
        )
        candidate_id = (
            f"candidate_{rank:02d}_"
            f"{int(row['alignment_start'])}-{int(row['alignment_end'])}"
        )
        candidate_dir = output_dir / "top_candidates" / candidate_id
        distance_matrix_for_candidate(row, discovery, context).to_csv(
            candidate_dir / f"{candidate_id}_pairwise_p_distance.csv"
        )
        export_diagnostic_positions(
            row,
            encoded,
            discovery,
            candidate_dir,
            candidate_id,
        )

    write_summary(
        output_dir,
        args.alignment.resolve(),
        args.metadata.resolve(),
        discovery,
        candidates,
        shortlisted,
    )

    print(f"Candidates passing filters: {len(candidates)}")
    print(
        "Diagnostic sites in all three contrasts: "
        f"{int(candidates['all_contrasts_have_diagnostic_sites'].sum())}"
    )
    print(
        "All three distance margins positive:     "
        f"{int(candidates['all_margins_positive'].sum())} (supporting result)"
    )
    print(f"Top candidates exported:    {len(shortlisted)}")
    print(f"Results written to:         {output_dir}")

if __name__ == "__main__":
    main()

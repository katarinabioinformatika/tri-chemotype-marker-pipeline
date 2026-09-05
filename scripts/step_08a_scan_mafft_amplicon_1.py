#!/usr/bin/env python3

"""Find Amplicon 1 candidates separating T-2/HT-2 from Type B genomes."""

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
    / "amplicon_1_t2_vs_typeb_v2_terminal_site_filter"
)

TARGET_LENGTHS_BP = (1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
MIN_INDIVIDUAL_LENGTH_FRACTION = 0.80
SCAN_STEP_COLUMNS = 100

# Search 300-column terminal zones for continuous 22-bp binding sites.
# Sites allow at most three mismatches per genome and none in the 3′-terminal 5 nt.
TERMINAL_SEARCH_COLUMNS = 300
BINDING_SITE_LENGTH = 22
MIN_BINDING_SITE_MEAN_CONSERVATION = 0.95
MAX_BINDING_SITE_MISMATCHES_PER_GENOME = 3
THREE_PRIME_BASES = 5
MAX_THREE_PRIME_MISMATCHES_PER_GENOME = 0

MIN_CANDIDATE_COVERAGE = 0.85

MIN_PAIR_COMPARABLE_FRACTION = 0.70
MIN_DIAGNOSTIC_GROUP_COVERAGE = 0.80
MIN_DIAGNOSTIC_CONSENSUS = 0.75
TOP_CANDIDATES_TO_EXPORT = 10

GENOME_ID_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)", re.IGNORECASE)
VALID_BASES = set("ACGT")

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find Amplicon 1 candidates for T-2/HT-2 versus Type B."
    )
    parser.add_argument("--alignment", type=Path, default=DEFAULT_ALIGNMENT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top", type=int, default=TOP_CANDIDATES_TO_EXPORT)
    return parser.parse_args()

def extract_genome_id(text: object) -> str | None:
    match = GENOME_ID_PATTERN.search(str(text))
    return match.group(1).upper() if match else None

def normalise_chemotype(value: object) -> str | None:
    """Map confirmed and asterisk-labelled records into four scan groups."""
    if pd.isna(value):
        return None
    text = str(value).strip().upper().replace("–", "-").replace("—", "-")
    compact = re.sub(r"[\s_*]", "", text)
    if compact in {"3-ADON", "3ADON"}:
        return "3-ADON"
    if compact in {"15-ADON", "15ADON"}:
        return "15-ADON"
    if compact == "NIV":
        return "NIV"
    if "T-2/HT-2" in compact or "T2/HT2" in compact:
        return "T-2/HT-2"
    return None

def safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip()).strip("_")

def load_alignment(path: Path) -> list[SeqRecord]:
    if not path.exists():
        sys.exit(f"ERROR: Alignment not found: {path}")
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        sys.exit(f"ERROR: No FASTA records found: {path}")
    if len({len(record.seq) for record in records}) != 1:
        sys.exit("ERROR: FASTA records do not all have the same aligned length.")

    genome_ids = [extract_genome_id(record.id) for record in records]
    if any(genome_id is None for genome_id in genome_ids):
        bad = [r.id for r, g in zip(records, genome_ids) if g is None]
        sys.exit(f"ERROR: Could not recover Genome_ID from: {bad}")
    duplicates = sorted({x for x in genome_ids if genome_ids.count(x) > 1})
    if duplicates:
        sys.exit(f"ERROR: Duplicate Genome_ID values in alignment: {duplicates}")
    return records

def load_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: Metadata not found: {path}")
    metadata = pd.read_csv(path, encoding="utf-8-sig")
    required = {"Genome_ID", "Species", "Chemotype"}
    missing = sorted(required.difference(metadata.columns))
    if missing:
        sys.exit(f"ERROR: Missing metadata columns: {missing}")

    metadata = metadata.copy()
    metadata["Genome_ID"] = metadata["Genome_ID"].astype(str).str.strip().str.upper()
    if "Locus_class" in metadata.columns:
        core = metadata[metadata["Locus_class"].eq("core_TRI_cluster")].copy()
        if not core.empty:
            metadata = core

    conflicts = []
    for genome_id, subset in metadata.groupby("Genome_ID"):
        if len(subset[["Species", "Chemotype"]].drop_duplicates()) > 1:
            conflicts.append(genome_id)
    if conflicts:
        sys.exit(f"ERROR: Conflicting metadata rows for: {conflicts}")

    metadata = metadata.drop_duplicates("Genome_ID", keep="first")
    metadata["Separation_group"] = metadata["Chemotype"].map(normalise_chemotype)
    metadata["Evidence_status"] = np.where(
        metadata["Chemotype"].astype(str).str.contains(r"\*", regex=True),
        "associated",
        "confirmed",
    )
    return metadata

def join_alignment_metadata(
    records: list[SeqRecord], metadata: pd.DataFrame
) -> tuple[list[SeqRecord], pd.DataFrame]:
    metadata_by_id = metadata.set_index("Genome_ID", drop=False)
    joined_rows = []
    labelled_records = []

    for record in records:
        genome_id = extract_genome_id(record.id)
        if genome_id not in metadata_by_id.index:
            sys.exit(f"ERROR: No metadata match for aligned genome: {genome_id}")
        row = metadata_by_id.loc[genome_id].to_dict()
        if row["Separation_group"] is None:
            continue
        row["FASTA_header"] = record.description
        joined_rows.append(row)
        labelled_records.append(record)

    joined = pd.DataFrame(joined_rows).reset_index(drop=True)
    required_groups = {"3-ADON", "15-ADON", "NIV", "T-2/HT-2"}
    missing = sorted(required_groups.difference(set(joined["Separation_group"])))
    if missing:
        sys.exit(f"ERROR: Missing labelled groups: {missing}")
    return labelled_records, joined

def encode_alignment(records: list[SeqRecord]) -> np.ndarray:
    lookup = np.full(256, 4, dtype=np.uint8)
    for code, base in enumerate(b"ACGT"):
        lookup[base] = code
        lookup[ord(chr(base).lower())] = code
    matrix = np.empty((len(records), len(records[0].seq)), dtype=np.uint8)
    for row, record in enumerate(records):
        raw = np.frombuffer(str(record.seq).encode("ascii"), dtype=np.uint8)
        matrix[row] = lookup[raw]
    return matrix

def cumulative_counts(values: np.ndarray) -> np.ndarray:
    prefix = np.zeros((values.shape[0], values.shape[1] + 1), dtype=np.int32)
    np.cumsum(values, axis=1, dtype=np.int32, out=prefix[:, 1:])
    return prefix

def column_prefixes(encoded: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    counts = np.stack([(encoded == base).sum(axis=0) for base in range(4)])
    valid = counts.sum(axis=0)
    coverage = valid / encoded.shape[0]
    conservation = np.divide(
        counts.max(axis=0), valid, out=np.zeros(valid.shape, float), where=valid > 0
    )
    return (
        np.concatenate(([0.0], np.cumsum(coverage))),
        np.concatenate(([0.0], np.cumsum(conservation))),
    )

def interval_mean(prefix: np.ndarray, start: int, end: int) -> float:
    return float((prefix[end] - prefix[start]) / (end - start))

def find_conserved_terminal_site(
    encoded: np.ndarray,
    zone_start: int,
    zone_end: int,
    reverse_site: bool,
) -> dict | None:
    """
    Find a continuous 22-bp binding-site candidate shared across every genome.

    This is a feasibility consensus, not a final primer. For a forward primer,
    the 3-prime end is the right side of the site. For a reverse primer, the
    reverse-complemented primer's 3-prime end corresponds to the left side of
    the forward-strand binding site.
    """
    best = None
    bases = np.asarray(list("ACGT"))
    last_start = zone_end - BINDING_SITE_LENGTH

    for site_start in range(zone_start, last_start + 1):
        site_end = site_start + BINDING_SITE_LENGTH
        window = encoded[:, site_start:site_end]

        if not np.all(window < 4):
            continue

        counts = np.stack([(window == base).sum(axis=0) for base in range(4)])
        consensus_codes = counts.argmax(axis=0)
        conservation_by_position = counts.max(axis=0) / window.shape[0]
        mean_conservation = float(conservation_by_position.mean())
        if mean_conservation < MIN_BINDING_SITE_MEAN_CONSERVATION:
            continue

        mismatch_matrix = window != consensus_codes[np.newaxis, :]
        mismatches_per_genome = mismatch_matrix.sum(axis=1)
        three_prime_slice = (
            mismatch_matrix[:, :THREE_PRIME_BASES]
            if reverse_site
            else mismatch_matrix[:, -THREE_PRIME_BASES:]
        )
        three_prime_mismatches = three_prime_slice.sum(axis=1)

        if mismatches_per_genome.max() > MAX_BINDING_SITE_MISMATCHES_PER_GENOME:
            continue
        if three_prime_mismatches.max() > MAX_THREE_PRIME_MISMATCHES_PER_GENOME:
            continue

        consensus_forward = "".join(bases[consensus_codes])
        oligo = (
            str(Seq(consensus_forward).reverse_complement())
            if reverse_site
            else consensus_forward
        )
        exact_fraction = float(np.mean(mismatches_per_genome == 0))
        candidate = {
            "alignment_start": site_start + 1,
            "alignment_end": site_end,
            "suggested_oligo_5to3": oligo,
            "mean_conservation": mean_conservation,
            "exact_fraction": exact_fraction,
            "max_mismatches_any_genome": int(mismatches_per_genome.max()),
            "max_3prime_mismatches_any_genome": int(
                three_prime_mismatches.max()
            ),
        }
        score = (
            exact_fraction,
            mean_conservation,
            -int(mismatches_per_genome.max()),
        )
        if best is None or score > best["_score"]:
            candidate["_score"] = score
            best = candidate

    if best is not None:
        best.pop("_score")
    return best

def median_physical_length(base_prefix: np.ndarray, start: int, end: int) -> float:
    return float(np.median(base_prefix[:, end] - base_prefix[:, start]))

def find_end_for_target(
    base_prefix: np.ndarray, start: int, target_bp: int, alignment_length: int
) -> int | None:
    if median_physical_length(base_prefix, start, alignment_length) < target_bp:
        return None
    low, high = start + 1, alignment_length
    while low < high:
        middle = (low + high) // 2
        if median_physical_length(base_prefix, start, middle) >= target_bp:
            high = middle
        else:
            low = middle + 1
    return low

def make_pair_data(encoded: np.ndarray):
    pair_i, pair_j = np.triu_indices(encoded.shape[0], k=1)
    comparable = (encoded[pair_i] < 4) & (encoded[pair_j] < 4)
    mismatch = comparable & (encoded[pair_i] != encoded[pair_j])
    return pair_i, pair_j, cumulative_counts(comparable), cumulative_counts(mismatch)

def pair_mask(
    groups: np.ndarray,
    pair_i: np.ndarray,
    pair_j: np.ndarray,
    group_a: set[str],
    group_b: set[str] | None = None,
) -> np.ndarray:
    ai = np.isin(groups[pair_i], list(group_a))
    aj = np.isin(groups[pair_j], list(group_a))
    if group_b is None:
        return ai & aj
    bi = np.isin(groups[pair_i], list(group_b))
    bj = np.isin(groups[pair_j], list(group_b))
    return (ai & bj) | (bi & aj)

def distance_margin(
    distances: np.ndarray,
    within_a: np.ndarray,
    within_b: np.ndarray,
    between: np.ndarray,
) -> tuple[float, float, float, float]:
    a = distances[within_a & np.isfinite(distances)]
    b = distances[within_b & np.isfinite(distances)]
    across = distances[between & np.isfinite(distances)]
    if not len(a) or not len(b) or not len(across):
        return (float("nan"),) * 4
    max_a = float(a.max())
    max_b = float(b.max())
    min_between = float(across.min())
    return max_a, max_b, min_between, min_between - max(max_a, max_b)

def diagnostic_sites(
    encoded: np.ndarray,
    groups: np.ndarray,
    start: int,
    end: int,
    group_a: set[str],
    group_b: set[str],
    name_a: str,
    name_b: str,
) -> pd.DataFrame:
    seq_a = encoded[np.isin(groups, list(group_a)), start:end]
    seq_b = encoded[np.isin(groups, list(group_b)), start:end]
    counts_a = np.stack([(seq_a == base).sum(axis=0) for base in range(4)])
    counts_b = np.stack([(seq_b == base).sum(axis=0) for base in range(4)])
    valid_a, valid_b = counts_a.sum(axis=0), counts_b.sum(axis=0)
    code_a, code_b = counts_a.argmax(axis=0), counts_b.argmax(axis=0)
    cov_a, cov_b = valid_a / seq_a.shape[0], valid_b / seq_b.shape[0]
    cons_a = np.divide(
        counts_a.max(axis=0), valid_a, out=np.zeros(valid_a.shape, float), where=valid_a > 0
    )
    cons_b = np.divide(
        counts_b.max(axis=0), valid_b, out=np.zeros(valid_b.shape, float), where=valid_b > 0
    )
    keep = (
        (cov_a >= MIN_DIAGNOSTIC_GROUP_COVERAGE)
        & (cov_b >= MIN_DIAGNOSTIC_GROUP_COVERAGE)
        & (cons_a >= MIN_DIAGNOSTIC_CONSENSUS)
        & (cons_b >= MIN_DIAGNOSTIC_CONSENSUS)
        & (code_a != code_b)
    )
    offsets = np.flatnonzero(keep)
    bases = np.asarray(list("ACGT"))
    return pd.DataFrame(
        {
            "alignment_position": start + offsets + 1,
            "group_a": name_a,
            "group_b": name_b,
            "group_a_consensus": bases[code_a[offsets]],
            "group_b_consensus": bases[code_b[offsets]],
            "group_a_coverage": cov_a[offsets],
            "group_b_coverage": cov_b[offsets],
            "group_a_consensus_fraction": cons_a[offsets],
            "group_b_consensus_fraction": cons_b[offsets],
        }
    )

SCAN_NAME = "Amplicon 1: T-2/HT-2 versus Type B"
SCAN_DEFINITION = {
    "included_groups": {"T-2/HT-2", "3-ADON", "15-ADON", "NIV"},
    "contrasts": [
        (
            "t2_typeb",
            {"T-2/HT-2"},
            {"3-ADON", "15-ADON", "NIV"},
            "T-2/HT-2",
            "Type B",
        )
    ],
}

def scan_one_amplicon(
    records: list[SeqRecord], metadata: pd.DataFrame, definition: dict
) -> tuple[pd.DataFrame, dict, np.ndarray]:
    encoded = encode_alignment(records)
    groups = metadata["Separation_group"].to_numpy()
    alignment_length = encoded.shape[1]
    base_prefix = cumulative_counts(encoded < 4)
    coverage_prefix, conservation_prefix = column_prefixes(encoded)
    pair_i, pair_j, comparable_prefix, mismatch_prefix = make_pair_data(encoded)

    contrast_masks = {}
    for key, group_a, group_b, _, _ in definition["contrasts"]:
        contrast_masks[key] = (
            pair_mask(groups, pair_i, pair_j, group_a),
            pair_mask(groups, pair_i, pair_j, group_b),
            pair_mask(groups, pair_i, pair_j, group_a, group_b),
        )

    rows = []
    seen = set()
    for start in range(0, alignment_length, SCAN_STEP_COLUMNS):
        for target_bp in TARGET_LENGTHS_BP:
            end = find_end_for_target(base_prefix, start, target_bp, alignment_length)
            if end is None or (start, end) in seen:
                continue
            seen.add((start, end))

            lengths = base_prefix[:, end] - base_prefix[:, start]
            median_length = float(np.median(lengths))
            if median_length < min(TARGET_LENGTHS_BP) or median_length > max(TARGET_LENGTHS_BP) + 25:
                continue
            if lengths.min() < median_length * MIN_INDIVIDUAL_LENGTH_FRACTION:
                continue

            full_coverage = interval_mean(coverage_prefix, start, end)
            if full_coverage < MIN_CANDIDATE_COVERAGE:
                continue
            left_end = min(start + TERMINAL_SEARCH_COLUMNS, end)
            right_start = max(start, end - TERMINAL_SEARCH_COLUMNS)
            left_cov = interval_mean(coverage_prefix, start, left_end)
            right_cov = interval_mean(coverage_prefix, right_start, end)
            left_cons = interval_mean(conservation_prefix, start, left_end)
            right_cons = interval_mean(conservation_prefix, right_start, end)

            left_site = find_conserved_terminal_site(
                encoded, start, left_end, reverse_site=False
            )
            right_site = find_conserved_terminal_site(
                encoded, right_start, end, reverse_site=True
            )
            if left_site is None or right_site is None:
                continue

            comparable = comparable_prefix[:, end] - comparable_prefix[:, start]
            mismatches = mismatch_prefix[:, end] - mismatch_prefix[:, start]
            shorter = np.minimum(lengths[pair_i], lengths[pair_j])
            required = np.maximum(
                1, np.ceil(shorter * MIN_PAIR_COMPARABLE_FRACTION).astype(int)
            )
            distances = np.divide(
                mismatches,
                comparable,
                out=np.full(comparable.shape, np.nan, dtype=float),
                where=comparable >= required,
            )

            result = {
                "alignment_start": start + 1,
                "alignment_end": end,
                "alignment_width_columns": end - start,
                "target_length_bp": target_bp,
                "median_length_bp": median_length,
                "min_length_bp": int(lengths.min()),
                "max_length_bp": int(lengths.max()),
                "q10_length_bp": float(np.quantile(lengths, 0.10)),
                "q90_length_bp": float(np.quantile(lengths, 0.90)),
                "mean_candidate_coverage": full_coverage,
                "left_end_coverage": left_cov,
                "right_end_coverage": right_cov,
                "left_end_conservation": left_cons,
                "right_end_conservation": right_cons,
                "left_site_start": left_site["alignment_start"],
                "left_site_end": left_site["alignment_end"],
                "left_site_consensus_5to3": left_site["suggested_oligo_5to3"],
                "left_site_mean_conservation": left_site["mean_conservation"],
                "left_site_exact_fraction": left_site["exact_fraction"],
                "left_site_max_mismatches": left_site["max_mismatches_any_genome"],
                "left_site_max_3prime_mismatches": left_site[
                    "max_3prime_mismatches_any_genome"
                ],
                "right_site_start": right_site["alignment_start"],
                "right_site_end": right_site["alignment_end"],
                "right_site_consensus_5to3": right_site["suggested_oligo_5to3"],
                "right_site_mean_conservation": right_site["mean_conservation"],
                "right_site_exact_fraction": right_site["exact_fraction"],
                "right_site_max_mismatches": right_site["max_mismatches_any_genome"],
                "right_site_max_3prime_mismatches": right_site[
                    "max_3prime_mismatches_any_genome"
                ],
                "minimum_pair_comparable_bp": int(comparable.min()),
                "mean_pair_comparable_bp": float(comparable.mean()),
            }

            site_counts = []
            margins = []
            for key, group_a, group_b, name_a, name_b in definition["contrasts"]:
                sites = diagnostic_sites(
                    encoded, groups, start, end, group_a, group_b, name_a, name_b
                )
                within_a, within_b, between = contrast_masks[key]
                margin_values = distance_margin(
                    distances, within_a, within_b, between
                )
                site_count = len(sites)
                site_counts.append(site_count)
                margins.append(margin_values[3])
                result[f"diagnostic_sites_{key}"] = site_count
                result[f"diagnostic_density_{key}_per_kb"] = (
                    site_count / (median_length / 1000.0)
                )
                result[f"within_{key}_group_a_max"] = margin_values[0]
                result[f"within_{key}_group_b_max"] = margin_values[1]
                result[f"between_{key}_min"] = margin_values[2]
                result[f"margin_{key}"] = margin_values[3]

            result["balanced_diagnostic_sites"] = min(site_counts)
            result["total_diagnostic_sites"] = sum(site_counts)
            result["balanced_diagnostic_density_per_kb"] = min(
                result[f"diagnostic_density_{key}_per_kb"]
                for key, *_ in definition["contrasts"]
            )
            finite_margins = np.asarray(margins, dtype=float)
            result["positive_margin_count"] = int(np.sum(finite_margins > 0))
            result["balanced_margin"] = (
                float(np.min(finite_margins))
                if np.all(np.isfinite(finite_margins))
                else np.nan
            )
            result["all_required_margins_positive"] = bool(
                np.all(finite_margins > 0)
            )
            rows.append(result)

    candidates = pd.DataFrame(rows)
    context = {
        "pair_i": pair_i,
        "pair_j": pair_j,
        "comparable_prefix": comparable_prefix,
        "mismatch_prefix": mismatch_prefix,
    }
    return candidates, context, encoded

def remove_overlapping(candidates: pd.DataFrame) -> pd.DataFrame:
    selected = []
    for index, row in candidates.iterrows():
        start, end = int(row["alignment_start"]), int(row["alignment_end"])
        width = end - start + 1
        if any(
            max(0, min(end, e) - max(start, s) + 1) / min(width, w) >= 0.80
            for s, e, w, _ in selected
        ):
            continue
        selected.append((start, end, width, index))
    return candidates.loc[[item[3] for item in selected]].copy()

def rank_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Prioritise primer compatibility, then balanced diagnostic signal."""
    if candidates.empty:
        return candidates
    ranked = candidates.sort_values(
        by=[
            "left_site_exact_fraction",
            "right_site_exact_fraction",
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
        ascending=[False, False, False, False, False, False, False, False, False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    ranked.insert(0, "overall_rank", np.arange(1, len(ranked) + 1))
    return ranked

def export_candidate(
    row: pd.Series,
    shortlist_rank: int,
    records: list[SeqRecord],
    metadata: pd.DataFrame,
    encoded: np.ndarray,
    definition: dict,
    output_dir: Path,
) -> None:
    start = int(row["alignment_start"]) - 1
    end = int(row["alignment_end"])
    candidate_id = f"candidate_{shortlist_rank:02d}_{start + 1}-{end}"
    candidate_dir = output_dir / "top_candidates" / candidate_id
    candidate_dir.mkdir(parents=True, exist_ok=True)

    aligned, ungapped = [], []
    for record, (_, meta) in zip(records, metadata.iterrows()):
        sequence = str(record.seq[start:end]).upper()
        sequence_ungapped = "".join(base for base in sequence if base in VALID_BASES)
        header = "|".join(
            [
                safe_name(meta["Separation_group"]),
                safe_name(meta["Species"]),
                meta["Genome_ID"],
                safe_name(meta["Evidence_status"]),
            ]
        )
        aligned.append(SeqRecord(Seq(sequence), id=header, description=""))
        ungapped.append(SeqRecord(Seq(sequence_ungapped), id=header, description=""))

    SeqIO.write(aligned, candidate_dir / f"{candidate_id}_aligned.fasta", "fasta")
    SeqIO.write(
        ungapped,
        candidate_dir / f"{candidate_id}_ungapped_for_Geneious.fasta",
        "fasta",
    )

    groups = metadata["Separation_group"].to_numpy()
    tables = []
    for key, group_a, group_b, name_a, name_b in definition["contrasts"]:
        table = diagnostic_sites(
            encoded, groups, start, end, group_a, group_b, name_a, name_b
        )
        table.insert(0, "contrast", key)
        tables.append(table)
    pd.concat(tables, ignore_index=True).to_csv(
        candidate_dir / f"{candidate_id}_diagnostic_positions.csv", index=False
    )

def write_scan_summary(
    scan_name: str,
    metadata: pd.DataFrame,
    ranked: pd.DataFrame,
    shortlisted: pd.DataFrame,
    output_dir: Path,
) -> None:
    lines = [
        scan_name,
        "=" * len(scan_name),
        "",
        f"Genomes used: {len(metadata)}",
        "",
        "Groups:",
    ]
    for group, count in metadata["Separation_group"].value_counts().sort_index().items():
        lines.append(f"  {group:12s}: {count}")
    lines.extend(
        [
            "",
            f"Candidates passing filters: {len(ranked)}",
            f"Candidates with all required margins positive: "
            f"{int(ranked['all_required_margins_positive'].sum()) if len(ranked) else 0}",
            f"Non-redundant candidates exported: {len(shortlisted)}",
            "",
            "The terminal regions are conserved search zones, not designed primers.",
            "Primer design and cross-genome testing must be completed in Geneious Prime.",
        ]
    )
    (output_dir / "analysis_summary.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

def main() -> None:
    args = parse_arguments()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    all_records = load_alignment(args.alignment.resolve())
    metadata = load_metadata(args.metadata.resolve())
    labelled_records, labelled_metadata = join_alignment_metadata(
        all_records, metadata
    )
    record_by_id = {
        extract_genome_id(record.id): record for record in labelled_records
    }

    scan_metadata = labelled_metadata[
        labelled_metadata["Separation_group"].isin(
            SCAN_DEFINITION["included_groups"]
        )
    ].copy().reset_index(drop=True)
    scan_records = [record_by_id[x] for x in scan_metadata["Genome_ID"]]
    scan_metadata.to_csv(output_dir / "genomes_used.csv", index=False)

    print(SCAN_NAME)
    print("-" * len(SCAN_NAME))
    print(f"Alignment records: {len(all_records)}")
    print(f"Alignment length:  {len(all_records[0].seq):,} columns")
    print(f"Genomes used:      {len(scan_metadata)}")
    print()
    print(scan_metadata["Separation_group"].value_counts().sort_index().to_string())
    print("\nScanning candidate regions...")

    candidates, _, encoded = scan_one_amplicon(
        scan_records, scan_metadata, SCAN_DEFINITION
    )
    ranked = rank_candidates(candidates)
    ranked.to_csv(output_dir / "all_ranked_candidates.csv", index=False)

    if ranked.empty:
        shortlisted = ranked.copy()
        print("No candidates passed the current filters.")
    else:
        shortlisted = remove_overlapping(ranked).head(max(1, args.top)).copy()
        shortlisted = shortlisted.reset_index(drop=True)
        shortlisted.insert(
            0, "shortlist_rank", np.arange(1, len(shortlisted) + 1)
        )
        shortlisted.to_csv(
            output_dir / "top_nonredundant_candidates.csv", index=False
        )
        for _, row in shortlisted.iterrows():
            export_candidate(
                row,
                int(row["shortlist_rank"]),
                scan_records,
                scan_metadata,
                encoded,
                SCAN_DEFINITION,
                output_dir,
            )
        print(f"Candidates passing filters: {len(ranked)}")
        print(f"Top candidates exported:   {len(shortlisted)}")

    write_scan_summary(
        SCAN_NAME,
        scan_metadata,
        ranked,
        shortlisted,
        output_dir,
    )
    print(f"Results written to: {output_dir}")

if __name__ == "__main__":
    main()

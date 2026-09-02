"""
part3_db_conclusion.py — Part 3 of the Deletion Bridge Analysis (see
DELETION_BRIDGE_ANALYSIS_PLAN.md at the project root — read that file
first, it is the canonical, locked spec, decisions #5-6 in particular).

For every mmc5 chain, one row per deletion-bridge pair (or a single
placeholder row for chains with no DB pair at all), recording whether
the segments strictly inside that DB span are ALTERED in the final
structure (Part 2) relative to the original Baca chain (Part 1), and if
so, how.

Per-segment verdict (using Piece.moved / Piece.flipped, already computed
by final_structure_assembly.assemble_final_structure -- not re-derived
here):
  UNCHANGED        moved=False, flipped=False -- exact same position,
                    same neighbors, same orientation as the original chain.
  FLIPPED-IN-PLACE  moved=False, flipped=True  -- stayed adjacent to its
                    original native neighbors, only its own internal
                    orientation reversed (the apostrophe case, no relocation).
  MOVED             moved=True (regardless of flipped) -- relocated to a
                    genuinely different neighborhood; the per-segment
                    detail column separately notes if it's also flipped.
  DELETED           reserved for a segment genuinely unaccounted for
                    anywhere in the final structure. VERIFIED EMPIRICALLY
                    (2026-08-13) that this can never actually occur for
                    this dataset: n_unresolved_rearr_ends == 0 across all
                    194 DB-bearing chains (mmc5's own Rearrangement-number
                    groups never cross chain boundaries -- already
                    established project-wide), and the assembly algorithm
                    itself never destroys a piece, only relocates/
                    reorients it. Kept in the schema per the locked
                    decision (a chain with an unresolved rearrangement end
                    WOULD produce it, this cohort simply never triggers
                    that case) -- expected count is 0, and the run below
                    confirms it, not assumes it.

DB-pair-level (row) verdict: `Altered? (Yes/No)` is No only if EVERY
segment in that DB pair is UNCHANGED; otherwise Yes, with a `sub_type`
column reporting the set of distinct non-UNCHANGED verdicts found across
the pair's segments (comma-joined if more than one distinct type
occurs -- decision #6's "or a combination where applicable"). A
`per_segment_detail` column gives the full per-letter breakdown for
direct cross-verification against the Part 1 / Part 2 images, per the
plan's required cross-verification step.

Output: results/Deletion Bridge Analysis on Baca Chains/
  DB Conclusion/db_conclusion.csv
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from db_segments_common import (  # noqa: E402
    load_mmc5,
    get_all_db_pairs_with_segments,
    chain_has_db,
)
from final_structure_assembly import (  # noqa: E402
    assemble_final_structure,
    build_piece_lookup,
)

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
TOP_FOLDER = os.path.join(RESULTS_ROOT, "Deletion Bridge Analysis on Baca Chains")
PART3_FOLDER = os.path.join(TOP_FOLDER, "DB Conclusion")
OUTPUT_CSV = os.path.join(PART3_FOLDER, "db_conclusion.csv")


def _segment_verdict(piece):
    if piece.moved:
        return 'MOVED'
    if piece.flipped:
        return 'FLIPPED-IN-PLACE'
    return 'UNCHANGED'


def compute_db_pair_rows(patient_id, chain_number, chain_rows):
    """Returns a list of row dicts, one per DB pair in this chain (empty
    list if the chain has no DB pair -- the caller emits the single
    N-row placeholder in that case)."""
    db_pairs = get_all_db_pairs_with_segments(chain_rows)
    if not db_pairs:
        return []

    result = assemble_final_structure(chain_rows)
    piece_lookup = build_piece_lookup(result)

    rows = []
    for dbp in db_pairs:
        chrom = dbp['chromosome']
        per_segment = []
        verdicts = []
        for seg in dbp['segments']:
            key = (chrom, seg['low_bp'], seg['high_bp'])
            piece = piece_lookup[key]  # must exist -- see module docstring
            verdict = _segment_verdict(piece)
            verdicts.append(verdict)
            flip_note = ' (flipped)' if (piece.moved and piece.flipped) else ''
            per_segment.append(f"{seg['letter']}:{verdict}{flip_note}")

        distinct_non_unchanged = sorted(set(v for v in verdicts if v != 'UNCHANGED'))
        altered = 'Yes' if distinct_non_unchanged else 'No'
        sub_type = ', '.join(distinct_non_unchanged) if distinct_non_unchanged else ''

        rows.append({
            'patient_id': patient_id,
            'chain_number': chain_number,
            'db_exists_in_chain': 'Y',
            'db_index': dbp['db_index'],
            'db_anchor_low_bp': dbp['low_bp'],
            'db_anchor_high_bp': dbp['high_bp'],
            'chromosome': chrom,
            'segment_letters_covered': ','.join(s['letter'] for s in dbp['segments']),
            'n_segments': len(dbp['segments']),
            'altered': altered,
            'sub_type': sub_type,
            'per_segment_detail': '; '.join(per_segment),
        })
    return rows


CSV_COLUMNS = [
    'patient_id', 'chain_number', 'db_exists_in_chain', 'db_index',
    'db_anchor_low_bp', 'db_anchor_high_bp', 'chromosome',
    'segment_letters_covered', 'n_segments', 'altered', 'sub_type',
    'per_segment_detail',
]


def run_all(output_csv=OUTPUT_CSV):
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df = load_mmc5()

    all_rows = []
    n_chains, n_db_rows, n_no_db_rows, n_failed = 0, 0, 0, 0
    for (patient_id, chain_number), grp in df.groupby(['Individual', 'Chain number']):
        n_chains += 1
        try:
            if not chain_has_db(grp):
                all_rows.append({
                    'patient_id': patient_id, 'chain_number': chain_number,
                    'db_exists_in_chain': 'N',
                    'db_index': '', 'db_anchor_low_bp': '', 'db_anchor_high_bp': '',
                    'chromosome': '', 'segment_letters_covered': '', 'n_segments': '',
                    'altered': '', 'sub_type': '', 'per_segment_detail': '',
                })
                n_no_db_rows += 1
                continue
            rows = compute_db_pair_rows(patient_id, chain_number, grp)
            all_rows.extend(rows)
            n_db_rows += len(rows)
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    n_deleted = sum(1 for r in all_rows if r['sub_type'] and 'DELETED' in r['sub_type'])
    print(f'Chains processed: {n_chains}  Failed: {n_failed}')
    print(f'DB-pair rows (Y): {n_db_rows}  No-DB placeholder rows (N): {n_no_db_rows}')
    print(f'Total rows written: {len(all_rows)}  -> {output_csv}')
    print(f'DELETED occurrences: {n_deleted} (expected 0 -- see module docstring)')


if __name__ == '__main__':
    run_all()

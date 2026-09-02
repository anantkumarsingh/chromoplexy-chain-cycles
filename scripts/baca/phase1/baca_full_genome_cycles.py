"""
baca_full_genome_cycles.py  —  Phase 1: full-genome AMG cycle completion

Professor Arsuaga's 2026-06-18 directive: "complete the cycles" for the
Baca patients before adding biological constraints, so we can see what
extra information that provides and identify which patients actually have
chromoplexy.

compute_erg_cycle_structures() in baca_aberration_clinical_analysis.py only
ever sees the ERG-anchored chain (_get_erg_chain_rows: anchor rows + one-hop
inter_chr expansion). Most "open paths" reported there are an artifact of
that restriction — the dangling end's true rearrangement partner is usually
recorded elsewhere in the same patient's full CSV. Completing the cycles
means building the AMG from EVERY rearrangement row for the patient, genome
wide, for all 57 patients (not just the 26 ERG+ ones used by the chain
analysis) — this is what actually answers "does this patient have
chromoplexy", rather than the k_local heuristic used in
get_erg_chain_details_v2.

This module reuses the low-level AMG primitives from
baca_aberration_clinical_analysis.py (_build_amg_nodes_and_edges) and the
drawing primitives from baca_aberration_cycle_drawing.py (draw_full_chain)
rather than duplicating them — only the genome-wide-specific logic lives
here.

Chromoplexy flag definition: a patient is flagged if at least one CLOSED
cycle spans >= min_chromosomes (default 3) distinct chromosomes — i.e. a
single coordinated multi-chromosome rearrangement, which is the actual
chromoplexy phenomenon. A cycle spanning only 2 chromosomes is just a
simple reciprocal translocation, not chromoplexy, so min_chromosomes=2
over-counts. This is deliberately stricter than (and not directly
comparable to) Baca's ERG-specific 15/26 — that number answers "is the ERG
fusion embedded in a chromoplexy chain" among ERG+ patients only; this flag
answers "does ANY part of this patient's genome show a closed, coordinated,
>=3-chromosome cycle", across all 57 patients.

Outputs:
    results/full_genome_cycle_summary.csv          — one row per patient, full cycle structure string
    results/full_genome_cycle_compact_summary.csv  — one row per patient, compact stats for analysis
    results/full_genome_cycle_diagrams/*.png        — one single, expanded full-chain panel per patient
"""

import math
import os
import sys
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
from baca_aberration_clinical_analysis import (
    load_data,
    _build_amg_nodes_and_edges,
    BACA_CHROM_ABER_FILE_PATH,
    BACA_CLINICAL_PHENO_FILE_PATH,
)
from baca_aberration_cycle_drawing import (
    build_norm_map,
    node_x,
    draw_chr_track,
    _arc,
    _assign_rads,
    COL,
    CHR_SEP,
)

RESULTS_DIR = (
    '/Users/anantkumarsingh/projects/prostate_cancer'
    '/nih-tcga-prad/results'
)
FULL_GENOME_RESULTS_DIR = os.path.join(RESULTS_DIR, 'full_genome_chain_diagrams')
FULL_GENOME_SUMMARY_CSV = os.path.join(RESULTS_DIR, 'full_genome_cycle_summary.csv')
FULL_GENOME_COMPACT_SUMMARY_CSV = os.path.join(RESULTS_DIR, 'full_genome_cycle_compact_summary.csv')


# ── Data scoping ─────────────────────────────────────────────────────────────

def _get_full_patient_rows(patient_id, chrom_df):
    """
    Return every valid rearrangement row for this patient — the complete
    genome-wide rearrangement set, not just the ERG-anchored chain.
    Drops the small number of rows with missing breakpoint data (9/5719
    rows in the current dataset — a few fully-empty trailing CSV rows).
    """
    patient = chrom_df[chrom_df['Individual'] == patient_id].copy()
    required = [
        'Breakpoint 1 chromosome', 'Breakpoint 1 position', 'Breakpoint 1 strand',
        'Breakpoint 2 chromosome', 'Breakpoint 2 position', 'Breakpoint 2 strand',
    ]
    return patient.dropna(subset=required)


# ── Cycle traversal with chromosome tracking ────────────────────────────────

def _find_amg_cycles_with_chroms(all_nodes, rearr_adj, ref_adj):
    """
    Same alternating-edge traversal as _find_amg_cycles, but also records
    the set of chromosomes visited by each closed cycle. Needed to
    classify chromoplexy directly from the real cycle structure (a cycle
    spanning many chromosomes = chromoplexy) instead of the k_local
    heuristic.

    Returns:
        cycles     : list of (cycle_length, set_of_chromosomes)
        open_paths : list of lists (unchanged semantics vs _find_amg_cycles)
    """
    visited = set()
    cycles = []
    open_paths = []

    for start_node in sorted(all_nodes):
        if start_node in visited:
            continue
        if start_node not in rearr_adj:
            continue

        path = []
        current = start_node
        step = 'rearrangement'
        rearr_hops = 0
        is_open = False

        while current not in visited:
            visited.add(current)
            path.append(current)

            if step == 'rearrangement':
                nxt = rearr_adj.get(current)
                step = 'reference'
                if nxt is None:
                    is_open = True; break
                rearr_hops += 1
            else:
                nxt = ref_adj.get(current)
                step = 'rearrangement'
                if nxt is None or nxt == 'TELOMERE':
                    is_open = True; break

            current = nxt

        if (not is_open) and (current == start_node) and rearr_hops > 0:
            chroms = {n[0] for n in path}
            cycles.append((len(path), chroms))
        else:
            open_paths.append(path)

    return cycles, open_paths


def classify_chromoplexy_from_cycles(cycles_with_chroms, min_chromosomes=3):
    """
    A patient's full-genome cycle structure shows chromoplexy if at least
    one CLOSED cycle spans >= min_chromosomes distinct chromosomes — i.e.
    a single coordinated multi-chromosome rearrangement event, not just a
    local intra-chromosomal C2. This replaces the k_local>2-AND-chain_size>1
    heuristic with a direct readout of the actual cycle structure.

    min_chromosomes=3 is a first-pass, conservative threshold — only 4/57
    patients meet it. Worth recalibrating once we compare against Baca's
    own chromoplexy calls; min_chromosomes=2 (any closed cycle that crosses
    chromosomes at all) flags 35/57 and may be closer to their definition.
    """
    max_span = max((len(chroms) for _, chroms in cycles_with_chroms), default=0)
    return {
        'has_chromoplexy':      max_span >= min_chromosomes,
        'max_cycle_chrom_span': max_span,
        'n_closed_cycles':      len(cycles_with_chroms),
    }


# ── Per-patient computation ──────────────────────────────────────────────────

def compute_full_genome_cycle_structure(patient_id, chrom_df):
    """
    Build the complete AMG from ALL of a patient's rearrangement rows
    (genome-wide) and find its cycle structure. No counterfactual
    enumeration here — with ~100 rejoins per patient on average,
    (2n-1)!! is intractable; enumeration stays scoped to the ERG-chain
    analysis in compute_erg_cycle_structures().
    """
    patient_rows = _get_full_patient_rows(patient_id, chrom_df)
    if len(patient_rows) == 0:
        return None

    n_rejoins = len(patient_rows)
    n_dsbs = n_rejoins * 2

    all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(patient_rows)

    b_i = {chrom: len(positions) for chrom, positions in chr_to_dsb_pos.items()}
    k = len(b_i)
    max_b_i = max(b_i.values()) if b_i else 0

    cycles_with_chroms, open_paths = _find_amg_cycles_with_chroms(all_nodes, rearr_adj, ref_adj)
    cycle_lengths = sorted((c[0] for c in cycles_with_chroms), reverse=True)
    observed_str = ' + '.join(f'C{m}' for m in cycle_lengths) or 'none (all open paths)'
    observed_rank = (n_dsbs - len(cycle_lengths)) if cycle_lengths else None

    chromoplexy = classify_chromoplexy_from_cycles(cycles_with_chroms)

    return {
        'patient_id':            patient_id,
        'k':                     k,
        'n_rejoins':             n_rejoins,
        'n_dsbs':                n_dsbs,
        'max_b_i':               max_b_i,
        'observed_structure':    observed_str,
        'observed_n_cycles':     len(cycle_lengths),
        'observed_rank':         observed_rank,
        'open_paths':            len(open_paths),
        'has_chromoplexy':       chromoplexy['has_chromoplexy'],
        'max_cycle_chrom_span':  chromoplexy['max_cycle_chrom_span'],
    }


def run_full_genome_cycle_analysis(chrom_df, save_path=None):
    """
    Run compute_full_genome_cycle_structure for every patient in chrom_df.
    Optionally saves a summary CSV. Returns the summary DataFrame.
    """
    results = []
    for pid in sorted(chrom_df['Individual'].dropna().unique(), key=str):
        res = compute_full_genome_cycle_structure(pid, chrom_df)
        if res is not None:
            results.append(res)

    summary_df = pd.DataFrame(results)
    if save_path:
        summary_df.to_csv(save_path, index=False)
        print(f"Saved full-genome cycle summary -> {save_path} ({len(summary_df)} patients)")
    return summary_df


# ── Compact per-patient summary ──────────────────────────────────────────────
#
# observed_structure in compute_full_genome_cycle_structure can be hundreds
# of terms long for large patients (e.g. 240 cycles for PR-08-556) — fine
# for downstream processing, unreadable for review. This gives every field
# you'd otherwise have to read off the single-panel diagram, in one row.

def summarize_full_genome_cycles(patient_id, chrom_df, min_chromosomes=3):
    """
    Compact per-patient summary: everything you'd read off the full-chain
    diagram, condensed into one row.
    """
    patient_rows = _get_full_patient_rows(patient_id, chrom_df)
    if len(patient_rows) == 0:
        return None

    n_rejoins = len(patient_rows)
    n_dsbs = n_rejoins * 2

    all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(patient_rows)
    b_i = {chrom: len(positions) for chrom, positions in chr_to_dsb_pos.items()}
    k = len(b_i)
    max_b_i = max(b_i.values()) if b_i else 0

    cycles_with_chroms, open_paths = _find_amg_cycles_with_chroms(all_nodes, rearr_adj, ref_adj)
    cycle_lengths = [c[0] for c in cycles_with_chroms]
    chromoplexy = classify_chromoplexy_from_cycles(cycles_with_chroms, min_chromosomes)

    # Cycle-size histogram, e.g. "C2x30 + C4x10 + C6x5" — same information
    # as the raw structure string, far more readable for large patients.
    size_counts = Counter(cycle_lengths)
    histogram = ' + '.join(
        f'C{size}x{cnt}' for size, cnt in sorted(size_counts.items(), reverse=True)
    ) or 'none'

    # Fraction of this patient's rearrangement edges that ended up inside a
    # closed cycle vs. left dangling in an open path — a direct measure of
    # how much of the genome-wide AMG "completing the cycles" actually closes.
    rearr_edges_in_cycles = sum(cycle_lengths) // 2
    frac_rejoins_in_cycles = rearr_edges_in_cycles / n_rejoins if n_rejoins else 0.0

    pct_inter = (
        100.0 * (patient_rows['Class'] == 'inter_chr').sum() / n_rejoins
        if n_rejoins else 0.0
    )

    rank = (n_dsbs - len(cycle_lengths)) if cycle_lengths else None
    largest_cycle = max(cycle_lengths) if cycle_lengths else 0

    return {
        'patient_id':              patient_id,
        'k':                       k,
        'n_rejoins':               n_rejoins,
        'n_dsbs':                  n_dsbs,
        'max_b_i':                 max_b_i,
        'pct_inter':               round(pct_inter, 1),
        'n_closed_cycles':         len(cycle_lengths),
        'open_paths':              len(open_paths),
        'frac_rejoins_in_cycles':  round(frac_rejoins_in_cycles, 3),
        'rank':                    rank,
        'largest_cycle':           largest_cycle,
        'max_cycle_chrom_span':    chromoplexy['max_cycle_chrom_span'],
        'has_chromoplexy':         chromoplexy['has_chromoplexy'],
        'cycle_size_histogram':    histogram,
    }


def run_full_genome_compact_summary(chrom_df, save_path=None, min_chromosomes=3):
    """
    Run summarize_full_genome_cycles for every patient in chrom_df.
    Optionally saves a CSV. Returns the summary DataFrame.
    """
    results = []
    for pid in sorted(chrom_df['Individual'].dropna().unique(), key=str):
        res = summarize_full_genome_cycles(pid, chrom_df, min_chromosomes)
        if res is not None:
            results.append(res)

    summary_df = pd.DataFrame(results)
    if save_path:
        summary_df.to_csv(save_path, index=False)
        print(f"Saved compact cycle summary -> {save_path} ({len(summary_df)} patients)")
    return summary_df


# ── Visualization ─────────────────────────────────────────────────────────────
#
# A single, large Full Genome Chain panel per patient — all chromosomes, all
# rearrangement arcs, built from the patient's COMPLETE rearrangement set
# instead of the ERG-anchored chain. No clinical_df / ERG anchor needed —
# applies to all 57 patients.
#
# Every breakpoint gets a stable running number (BP1, BP2, ... sorted by
# chromosome then position) labeled directly on the diagram, plus a corner
# "BP # → Location" reference table listing each one's full genomic location
# + strand — the same feature added to the mmc5 Baca-chain diagrams in
# scripts/baca/mmc5/baca_chain_visualization.py (_bp_table_dims /
# _add_bp_location_table), ported here per user request. Duplicated locally
# rather than cross-imported from the mmc5 script, matching this project's
# established per-script self-containment convention (see CLAUDE_CONTEXT.md
# "scripts/baca/ layout"). IMPORTANT DIFFERENCE from mmc5: the raw
# chrom_aberrations_baca.csv data has no per-breakpoint ID of its own (mmc5
# has a real 'Breakpoint number' column from Baca's own ChainFinder output;
# chrom_aberrations_baca.csv only has a per-REARRANGEMENT 'Number' column) —
# so the BP numbering here is our OWN invented, deterministic scheme, not
# one of Baca's own IDs. This is called out in the on-diagram title too.
#
# This drawing function is self-contained (does not reuse the shared
# draw_full_chain/draw_bp_labels from core/baca_aberration_cycle_drawing.py,
# which are still used unmodified by the ERG-chain / mmc5 image families) so
# that adding the BP-numbering + table here cannot affect any other image
# set that shares those primitives.

_TABLE_Y_TOP, _TABLE_Y_BOTTOM = 0.87, 0.05  # figure-fraction vertical budget below the title


def _bp_table_dims(n, fig_h, max_col_groups=3):
    """
    Width/height (figure-fraction) the corner table needs for n entries.
    Same sizing logic as mmc5/baca_chain_visualization.py's _bp_table_dims:
    prefer growing the table DOWNWARD (one tall column) over sideways, only
    wrapping into extra columns when a single column can't fit at a
    readable font size — full-genome chains run even taller than mmc5
    chains for large patients (e.g. PR-08-556, 24 chromosomes), so the same
    "mostly vertical room" assumption holds even more strongly here.
    """
    fontsize = 8.0 if n <= 40 else (7.0 if n <= 80 else (6.0 if n <= 150 else 5.5))
    line_height_in = fontsize * 1.6 / 72.0
    title_in = 0.28

    height_budget_in = (_TABLE_Y_TOP - _TABLE_Y_BOTTOM) * fig_h
    rows_that_fit = max(6, int((height_budget_in - title_in) / line_height_in))

    n_col_groups = max(1, min(max_col_groups, math.ceil(n / rows_that_fit)))
    max_rows_per_col = math.ceil(n / n_col_groups)

    table_w = min(0.115 * n_col_groups + 0.045, 0.30)
    table_h = min((max_rows_per_col * line_height_in + title_in) / fig_h, _TABLE_Y_TOP - _TABLE_Y_BOTTOM)
    return table_w, table_h, n_col_groups, max_rows_per_col, fontsize


def _add_bp_location_table(fig, entries, table_w, table_h, n_col_groups, max_rows_per_col, fontsize):
    """
    Corner reference table: BP number -> full genomic location + strand,
    sorted by BP number (matches the on-diagram labels). entries is a list
    of (bp_num, chrom, pos, strand_str) already sorted by bp_num. Placed in
    FIGURE-fraction coordinates in the space the caller already reserved
    (via _bp_table_dims + shrinking the main axes) so it never overlaps the
    chain diagram itself.
    """
    n = len(entries)
    ax_table = fig.add_axes(
        [1.0 - table_w + 0.01, max(_TABLE_Y_TOP - table_h, 0.02), table_w - 0.015, table_h]
    )
    ax_table.axis('off')
    ax_table.set_title(f'BP # → Location  (n={n})', fontsize=9.5, fontweight='bold',
                        loc='left', pad=3)

    col_width = 1.0 / n_col_groups
    for gi in range(n_col_groups):
        chunk = entries[gi * max_rows_per_col: (gi + 1) * max_rows_per_col]
        lines = [f'BP{num:<4} chr{chrom}:{pos:,.0f} ({strand})' for num, chrom, pos, strand in chunk]
        ax_table.text(
            gi * col_width, 1.0, '\n'.join(lines),
            transform=ax_table.transAxes, ha='left', va='top',
            fontsize=fontsize, family='monospace', linespacing=1.6,
        )


def _assign_bp_numbers(chr_to_dsb_pos):
    """
    Our own stable breakpoint numbering (BP1, BP2, ...), sorted by
    chromosome then position. Not one of Baca's own IDs — see module note
    above. Returns {(chrom, pos): bp_num}.
    """
    all_bps = [(c, p) for c in sorted(chr_to_dsb_pos) for p in chr_to_dsb_pos[c]]
    return {cp: i + 1 for i, cp in enumerate(all_bps)}


def draw_bp_number_labels(ax, chrom, y, sorted_pos, nmap, bp_num_map, rearr_adj, max_positions=40):
    """
    On-diagram label 'BP{n} (+/-)' per breakpoint, staggered above/below in
    position order to avoid overlap — mirrors
    mmc5/baca_chain_visualization.py's draw_bp_number_labels. Skipped
    entirely past max_positions (genome-wide chromosomes can carry 100+
    DSBs for large patients — same threshold already used by the shared
    core draw_bp_labels for this reason); the corner BP#-Location table
    always lists every breakpoint regardless of on-diagram density.
    """
    if len(sorted_pos) > max_positions:
        return
    for i, pos in enumerate(sorted_pos):
        x = nmap[pos]
        strands = [s for s in ('+', '-') if (chrom, pos, s) in rearr_adj]
        strand_tag = '/'.join(strands) if strands else '?'
        bp_num = bp_num_map[(chrom, pos)]
        above = (i % 2 == 0)
        ly = y + 0.46 if above else y - 0.46
        va = 'bottom' if above else 'top'
        ax.text(
            x, ly, f'BP{bp_num} ({strand_tag})', ha='center', va=va,
            fontsize=8.5, fontweight='bold', color='#1A1A1A', zorder=6,
            bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='#999999', lw=0.7, alpha=0.9),
        )


def draw_patient_full_genome(patient_id, chrom_df, save_dir=FULL_GENOME_RESULTS_DIR):
    os.makedirs(save_dir, exist_ok=True)

    patient_rows = _get_full_patient_rows(patient_id, chrom_df)
    if patient_rows.empty:
        print(f'[{patient_id}] No rearrangement rows.'); return

    all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(patient_rows)
    cycles_with_chroms, open_paths = _find_amg_cycles_with_chroms(all_nodes, rearr_adj, ref_adj)
    cycles      = [c[0] for c in cycles_with_chroms]
    chromoplexy = classify_chromoplexy_from_cycles(cycles_with_chroms)

    n_rejoins = len(patient_rows)
    n_dsbs    = n_rejoins * 2
    k         = len(chr_to_dsb_pos)
    b_i       = {c: len(v) for c, v in chr_to_dsb_pos.items()}
    cycle_str = (' + '.join(f'C{n}' for n in sorted(cycles, reverse=True))
                 or 'none')
    theta     = (f"Θ({k}, "
                 f"({', '.join(str(b_i[c]) for c in sorted(chr_to_dsb_pos))}))")
    rank      = str(n_dsbs - len(cycles)) if cycles else 'N/A'

    sorted_chroms = sorted(chr_to_dsb_pos)
    chr_y  = {c: (len(sorted_chroms) - 1 - i) * CHR_SEP
              for i, c in enumerate(sorted_chroms)}
    nmaps  = {c: build_norm_map(chr_to_dsb_pos[c]) for c in sorted_chroms}
    y_span = (len(sorted_chroms) - 1) * CHR_SEP

    bp_num_map = _assign_bp_numbers(chr_to_dsb_pos)
    n_breakpoints = len(bp_num_map)
    table_entries = [
        (num, int(c), p, '/'.join(s for s in ('+', '-') if (c, p, s) in rearr_adj) or '?')
        for (c, p), num in sorted(bp_num_map.items(), key=lambda kv: kv[1])
    ]

    # Rearrangement arcs (same collection logic as the previous draw_full_chain call)
    seen, arc_data = set(), []
    for u, v in rearr_adj.items():
        key = tuple(sorted([u, v]))
        if key in seen:
            continue
        seen.add(key)
        xa = node_x(u[1], u[2], nmaps[u[0]])
        xb = node_x(v[1], v[2], nmaps[v[0]])
        ya, yb = chr_y[u[0]], chr_y[v[0]]
        arc_data.append((xa, ya, xb, yb, u[0] != v[0]))
    rads = _assign_rads(arc_data, [a[4] for a in arc_data])

    # Sizing — same physical-inches worst-case-arc-sagitta protection used in
    # mmc5/baca_chain_visualization.py, needed here too (and more so): full-
    # genome chains can carry far deeper/wider intrachromosomal arcs (up to
    # hundreds of DSBs on one chromosome for large patients) than any mmc5
    # chain (max 118 breakpoints total).
    fig_w = 24.0
    fig_h0 = max(7.0, k * 3.8 + 3.0)
    y_range0 = y_span * 1.05 + 4.0
    scale_y0 = fig_h0 / max(y_range0, 1e-6)
    scale_x = fig_w / 1.12

    top_y, bottom_y = max(chr_y.values()), min(chr_y.values())

    def _worst_sagitta_at(row_y):
        return max(
            (abs(rad) * abs(xb - xa) * scale_x
             for (xa, ya, xb, yb, is_inter), rad in zip(arc_data, rads)
             if ya == row_y and not is_inter),
            default=0.0,
        )

    extra_in_top = _worst_sagitta_at(top_y) * 1.15
    extra_in_bottom = _worst_sagitta_at(bottom_y) * 1.15
    extra_data_top = extra_in_top / scale_y0
    extra_data_bottom = extra_in_bottom / scale_y0
    fig_h = fig_h0 + extra_in_top + extra_in_bottom

    # Corner BP#-Location table dimensions, computed now that fig_h is
    # final, so the main diagram's right margin can be shrunk to make room
    # for it instead of the table overlapping the chain arcs.
    table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize = _bp_table_dims(n_breakpoints, fig_h)

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    x_margin = 0.06
    ax.set_xlim(-x_margin, 1.0 + x_margin)
    ax.set_ylim(bottom_y - 2.0 - extra_data_bottom, top_y + 2.0 + extra_data_top)
    ax.axis('off')

    chromoplexy_tag = 'CHROMOPLEXY' if chromoplexy['has_chromoplexy'] else 'no chromoplexy'
    ax.set_title(
        f'{patient_id} — Full Genome Chain\n'
        f'{theta}   n_breakpoints={n_breakpoints}   n_rearrangement_edges={n_rejoins}   '
        f'Observed: {cycle_str}  ({len(open_paths)} open)   Rank={rank}   '
        f'{chromoplexy_tag} (max cycle span={chromoplexy["max_cycle_chrom_span"]} chr)',
        fontsize=12, fontweight='bold', pad=14, loc='center',
    )

    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c])
        draw_bp_number_labels(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c], bp_num_map, rearr_adj)

    for (xa, ya, xb, yb, is_inter), rad in zip(arc_data, rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 1.8, '--', 0.85, '-', 4)

    # Legend at figure level (bottom-center), not ax.legend(loc=...) inside
    # the data axes — the established correct pattern for this drawing
    # family (see CLAUDE.md's chain-visualization margin bug-fix note).
    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Intrachromosomal rearrangement'),
        mpatches.Patch(color=COL['rearr_inter'], label='Interchromosomal rearrangement'),
    ], loc='lower center', ncol=2, fontsize=10, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    # Right margin shrunk by the table's reserved width so the corner
    # BP#-Location table sits in its own space, never overlapping the arcs.
    fig.subplots_adjust(left=0.04, right=max(0.96 - table_w, 0.6), top=0.90, bottom=0.08)
    _add_bp_location_table(fig, table_entries, table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize)

    out = os.path.join(save_dir, f'{patient_id}_full_genome.png')
    fig.savefig(out, dpi=150, facecolor='white')
    plt.close(fig)
    print(f'{patient_id} -> {out}  ({cycle_str})')


def run_all_patients_full_genome(chrom_df, save_dir=FULL_GENOME_RESULTS_DIR):
    """
    Run draw_patient_full_genome for every patient in chrom_df.
    """
    all_ids   = sorted(chrom_df['Individual'].dropna().unique(), key=str)
    n_total   = len(all_ids)
    succeeded = []
    failed    = []

    print(f"\n{'='*60}")
    print(f"  Full-Genome AMG Drawing  —  {n_total} patients")
    print(f"{'='*60}\n")

    for i, pid in enumerate(all_ids, 1):
        print(f"[{i:>2}/{n_total}] {pid}")
        try:
            draw_patient_full_genome(pid, chrom_df, save_dir=save_dir)
            succeeded.append(pid)
        except Exception as exc:
            print(f"         ERROR: {exc}")
            failed.append((pid, str(exc)))

    print(f"\n{'='*60}")
    print(f"  Done.  Succeeded: {len(succeeded)}/{n_total}")
    if failed:
        print(f"  Failed: {len(failed)}")
        for pid, err in failed:
            print(f"    {pid}: {err}")
    print(f"{'='*60}\n")
    return succeeded, failed


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    chrom_df, clinical_df = load_data(
        BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)
    run_full_genome_cycle_analysis(chrom_df, save_path=FULL_GENOME_SUMMARY_CSV)
    run_full_genome_compact_summary(chrom_df, save_path=FULL_GENOME_COMPACT_SUMMARY_CSV)
    run_all_patients_full_genome(chrom_df)

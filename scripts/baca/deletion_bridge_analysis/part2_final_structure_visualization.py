"""
part2_final_structure_visualization.py — Part 2 of the Deletion Bridge
Analysis (see DELETION_BRIDGE_ANALYSIS_PLAN.md at the project root —
read that file first, it is the canonical, locked spec).

REWRITTEN 2026-08-13 (three times): first to match Baca's own chain
style (one row per original chromosome, GAP markers for relocated
material — see final_structure_assembly.build_chromosome_rows), then
for a diamond/pairing/spacing polish pass, then a fourth-round polish
per further user feedback:
  - genuine solid-black TELOMERE bars are now drawn leading up to every
    outermost breakpoint (both truly isolated tips AND tips that got a
    real fusion, e.g. TELOMERE----BP2(-)) instead of a bare text label —
    see final_structure_assembly.build_chromosome_rows' `_render_tip`,
    which now always emits a ('telomere', node, has_real_fusion) row
    item first.
  - every real junction's diamond is now centered in a generously
    widened inter-cell gap (CELL_GAP raised substantially) instead of
    sitting almost on top of the adjacent breakpoint tick marks.
  - each chromosome's DB-span segment letters (A-Z) are drawn as their
    own clearly separated visual band, well below that row's own BP tick
    labels and well clear of the next chromosome row.
  - the deletion bridge itself is now drawn as a purple dotted arc
    connecting the DB anchor breakpoints' ACTUAL positions in the FINAL
    structure (same visual convention as Part 1's original-chain arc),
    so it's visible where the DB relationship's two ends physically
    ended up after rearrangement — tracked via a node-position map built
    while drawing each row.

Output: results/Deletion Bridge Analysis on Baca Chains/
  Final Rearrangement Structure/
    Final Structure After Rearrangements Images/
      {patient_id}_chain_{n}_final_rearranged_structure.png
      {patient_id}_chain_{n}_final_rearr_str_details.png
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'mmc5'))

from baca_chain_visualization import (  # noqa: E402  (reuse, don't duplicate)
    load_mmc5,
    DELETION_BRIDGE_COLOR,
)
from db_segments_common import (  # noqa: E402
    get_all_db_pairs_with_segments,
    chain_has_db,
)
from final_structure_assembly import (  # noqa: E402
    assemble_final_structure,
    classify_chain_rearrangements,
    build_chromosome_rows,
)
from part1_db_segments_visualization import DB_PALETTE, _db_color  # noqa: E402  (same colors, cross-image consistency)

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
TOP_FOLDER = os.path.join(RESULTS_ROOT, "Deletion Bridge Analysis on Baca Chains")
PART2_FOLDER = os.path.join(TOP_FOLDER, "Final Rearrangement Structure")
OUTPUT_DIR = os.path.join(PART2_FOLDER, "Final Structure After Rearrangements Images")

ROW_SEP = 11.0           # vertical distance between chromosome rows (raised for a dedicated segment-label band)
CELL_W = 0.20             # schematic width of one piece/gap/telomere cell
TIP_BAR_W = 0.14          # width of a telomere bar cell
CELL_GAP = 0.11           # generous breathing room between adjacent cells (diamonds live centered in here)
GAP_COLOR = '#EEEEEE'
GAP_HATCH_COLOR = '#BDBDBD'
TELOMERE_COLOR = '#111111'
COL_BP = '#C62828'
COL_JUNCTION = '#1565C0'

DB_LABEL_Y_OFFSET = -3.2   # segment-letter band, well clear of this row's own BP ticks AND the next row


def _build_db_label_lookup(db_pairs):
    lookup = {}
    use_prefix = len(db_pairs) > 1
    for dbp in db_pairs:
        for seg in dbp['segments']:
            key = (dbp['chromosome'], seg['low_bp'], seg['high_bp'])
            lookup[key] = (dbp['db_index'] if use_prefix else None, seg['letter'], _db_color(dbp['db_index']))
    return lookup


def _piece_label(piece, db_lookup):
    key = (piece.chrom, piece.low_bp, piece.high_bp)
    if key in db_lookup:
        db_index, letter, color = db_lookup[key]
        text = f'DB{db_index}:{letter}' if db_index is not None else letter
        if piece.flipped:
            text += "'"
        return text, color, True
    return None, '#9E9E9E', False


def _is_chain_terminal(end):
    return isinstance(end, tuple) and len(end) == 2 and end[0] == 'CHAIN_TERMINAL'


def _flatten_row_to_cells(row_items):
    """Expand a chromosome row (from build_chromosome_rows) into a flat
    list of draw cells. Every cell records whether a real fusion
    junction sits on its LEFT/RIGHT side (diamond needed) -- inter-piece
    junctions, and now also telomere-bar-to-piece junctions (a tip that
    got a real fusion)."""
    cells = []
    for item in row_items:
        kind = item[0]
        if kind == 'gap':
            cells.append({'kind': 'gap', 'meta': item[1]})
        elif kind == 'telomere':
            _, node, has_real_fusion = item
            cells.append({'kind': 'telomere', 'node': node, 'junction_right': has_real_fusion})
        else:
            _, mol, is_closed, reversed_disp = item
            pieces = list(reversed(mol.pieces)) if reversed_disp else list(mol.pieces)
            left_end = mol.right_free if reversed_disp else mol.left_free
            right_end = mol.left_free if reversed_disp else mol.right_free
            for i, piece in enumerate(pieces):
                is_first, is_last = (i == 0), (i == len(pieces) - 1)
                left_terminal = left_end if is_first else None
                right_terminal = right_end if is_last else None
                # Exactly ONE diamond per real junction: inter-piece
                # junctions are owned by the EARLIER piece's right side;
                # a piece's own left side only draws one when it's the
                # molecule's first piece connecting to a preceding
                # telomere-bar cell (drawn separately, not by this piece
                # -- see the 'telomere' cell's own junction_right above)
                # or a CHAIN_TERMINAL with no telomere-bar predecessor.
                cells.append({
                    'kind': 'piece', 'piece': piece, 'is_closed': is_closed,
                    'display_swapped': reversed_disp,
                    'junction_left': False,  # owned by the preceding telomere/piece cell instead
                    'junction_right': (not is_last) or _is_chain_terminal(right_terminal),
                    'left_terminal': left_terminal,
                    'right_terminal': right_terminal,
                })
    return cells


def _terminal_text(end):
    if end is None or not _is_chain_terminal(end):
        return None
    bp, strand = end[1]
    return f'→ BP{bp}({strand})'


def draw_chromosome_row(ax, y, chrom, cells, db_lookup, x_start, node_positions):
    """Draws one chromosome's row. Records the drawn (x, y) of every
    breakpoint tick mark into node_positions[bp] (list of x's, one per
    strand-side actually rendered) so the deletion-bridge arc can be
    drawn afterward, anchored to real, final positions."""
    x = x_start
    pending_junction = False  # True if the PREVIOUS cell wants a diamond drawn in the upcoming gap

    def _record(bp, xx):
        node_positions.setdefault(bp, []).append((xx, y))

    for cell in cells:
        if pending_junction:
            # `x` has already been advanced PAST the gap (to the next
            # cell's own x_lo) by the previous cell's own trailing
            # "x = x_hi + CELL_GAP" step, so the gap this diamond belongs
            # in is the interval [x - CELL_GAP, x], not [x, x + CELL_GAP]
            # -- bug found from the rendered image (2026-08-13): diamonds
            # were landing CELL_GAP too far right, inside the NEXT piece's
            # own gray bar instead of centered in the white space before
            # it.
            gap_mid = x - CELL_GAP / 2
            ax.plot(gap_mid, y, marker='D', markersize=11, color=COL_JUNCTION, zorder=8,
                    markeredgecolor='white', markeredgewidth=1.1)
            pending_junction = False

        if cell['kind'] == 'telomere':
            bp, strand = cell['node']
            x_lo, x_hi = x, x + TIP_BAR_W
            ax.plot([x_lo, x_hi], [y, y], color=TELOMERE_COLOR, linewidth=13, solid_capstyle='round', zorder=1)
            ax.plot([x_lo, x_lo], [y - 0.3, y + 0.3], color=TELOMERE_COLOR, linewidth=3, zorder=2)  # telomere end-cap
            ax.plot([x_hi, x_hi], [y - 0.26, y + 0.26], color=COL_BP, linewidth=2.6, zorder=3)
            ax.text((x_lo + x_hi) / 2, y + 0.9, 'TELOMERE', ha='center', va='bottom', fontsize=8,
                    color='#616161', style='italic')
            ax.text(x_hi, y - 0.55, f'BP{bp} ({strand})', ha='center', va='top', fontsize=9.5, fontweight='bold',
                    color='#1A1A1A', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#999999', lw=0.7, alpha=0.92))
            _record(bp, x_hi)
            x = x_hi + CELL_GAP
            pending_junction = cell['junction_right']
            continue

        if cell['kind'] == 'gap':
            x_lo, x_hi = x, x + CELL_W
            ax.add_patch(plt.Rectangle((x_lo, y - 0.38), x_hi - x_lo, 0.76, facecolor=GAP_COLOR,
                                        edgecolor=GAP_HATCH_COLOR, hatch='////', linewidth=0.9, zorder=1))
            ax.text((x_lo + x_hi) / 2, y, 'GAP', ha='center', va='center', fontsize=10,
                    color='#757575', fontweight='bold', zorder=2)
            meta = cell['meta']
            ax.text((x_lo + x_hi) / 2, y - 0.68, f"(was BP{meta['low_bp']}-BP{meta['high_bp']})",
                    ha='center', va='top', fontsize=7.5, color='#9E9E9E')
            x = x_hi + CELL_GAP
            continue

        # piece cell
        x_lo, x_hi = x, x + CELL_W
        piece = cell['piece']
        bar_color = '#B39DDB' if cell['is_closed'] else '#BDBDBD'
        pad = (x_hi - x_lo) * 0.08
        ax.plot([x_lo + pad, x_hi - pad], [y, y], color=bar_color, linewidth=13,
                solid_capstyle='butt', zorder=1)

        left_node, right_node = piece.display_nodes(cell['display_swapped'])
        for xx, node, side in [(x_lo + pad, left_node, 'L'), (x_hi - pad, right_node, 'R')]:
            ax.plot([xx, xx], [y - 0.26, y + 0.26], color=COL_BP, linewidth=2.6, zorder=3)
            bp, strand = node
            above = (side == 'L')
            ly = y + 0.55 if above else y - 0.55
            va = 'bottom' if above else 'top'
            ax.text(xx, ly, f'BP{bp} ({strand})', ha='center', va=va, fontsize=9.5, fontweight='bold',
                    color='#1A1A1A', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#999999', lw=0.7, alpha=0.92))
            _record(bp, xx)

        label, color, is_db = _piece_label(piece, db_lookup)
        mid = (x_lo + x_hi) / 2
        if is_db:
            ax.plot([x_lo + pad, x_hi - pad], [y + DB_LABEL_Y_OFFSET + 0.35, y + DB_LABEL_Y_OFFSET + 0.35],
                    color=color, linewidth=6, solid_capstyle='round', alpha=0.9, zorder=6)
            # connector line from the piece down to its own segment-label band, so it stays
            # legibly ATTACHED to the right piece even though the band is now far below the BP ticks
            ax.plot([mid, mid], [y - 0.65, y + DB_LABEL_Y_OFFSET + 0.35], color=color, linewidth=1.1,
                    alpha=0.5, zorder=5, linestyle=':')
            ax.text(mid, y + DB_LABEL_Y_OFFSET, label, ha='center', va='center', fontsize=11, fontweight='bold',
                    color=color, zorder=7,
                    bbox=dict(boxstyle='round,pad=0.22', fc='white', ec=color, lw=1.2, alpha=0.98))
        else:
            ax.text(mid, y + 0.02, f'BP{piece.low_bp}-BP{piece.high_bp}', ha='center', va='center',
                    fontsize=8, color='#BDBDBD', zorder=2)

        if cell['junction_right']:
            pending_junction = True

        x = x_hi + CELL_GAP

    return x  # final x extent used


def draw_deletion_bridge_arcs(ax, db_pairs, node_positions):
    """Purple dotted arc connecting each DB pair's two anchor breakpoints
    at their ACTUAL positions in the final structure (same color/style
    convention as Part 1's original-chain deletion-bridge arc) -- added
    2026-08-13 per user feedback: the final structure previously showed
    no indication of where the deletion bridge's two ends physically
    ended up after rearrangement. If an anchor breakpoint has more than
    one rendered strand-side (rare), the first recorded position is
    used -- sufficient to show WHERE the breakpoint physically sits."""
    for dbp in db_pairs:
        color = _db_color(dbp['db_index']) if len(db_pairs) > 1 else DELETION_BRIDGE_COLOR
        for bp in (dbp['low_bp'], dbp['high_bp']):
            if bp not in node_positions:
                return  # defensive: shouldn't happen, every breakpoint is drawn somewhere
        (xa, ya) = node_positions[dbp['low_bp']][0]
        (xb, yb) = node_positions[dbp['high_bp']][0]
        same_row = (ya == yb)
        rad = 0.3 if same_row else 0.15
        ax.add_patch(FancyArrowPatch(
            (xa, ya + 1.0 if same_row else ya), (xb, yb + 1.0 if same_row else yb),
            connectionstyle=f'arc3,rad={-rad}',
            arrowstyle='-', mutation_scale=15,
            color=DELETION_BRIDGE_COLOR, linewidth=2.2, linestyle=(0, (3, 2)), alpha=0.9, zorder=9,
        ))
        mid_x = (xa + xb) / 2
        mid_y = max(ya, yb) + (2.2 if same_row else 1.0)
        ax.text(mid_x, mid_y, f"DB{dbp['db_index']}: BP{dbp['low_bp']}↔BP{dbp['high_bp']}",
                ha='center', va='bottom', fontsize=9, fontweight='bold', color=DELETION_BRIDGE_COLOR,
                zorder=9, bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=DELETION_BRIDGE_COLOR, lw=0.8, alpha=0.92))


def draw_structure_diagram(patient_id, chain_number, chain_rows, db_pairs, db_lookup, result, save_dir):
    """Image 1 of the matched pair: the final structure diagram alone,
    generously spaced, no corner tables competing for room."""
    chrom_rows = build_chromosome_rows(chain_rows, result)
    sorted_chroms = sorted(chrom_rows.keys())
    cells_by_chrom = {c: _flatten_row_to_cells(chrom_rows[c]) for c in sorted_chroms}
    n_cells_max = max((len(c) for c in cells_by_chrom.values()), default=1)

    n_rows = len(sorted_chroms)
    fig_w = max(16.0, 2.5 + n_cells_max * (CELL_W + CELL_GAP) * 6.0)
    fig_h = max(7.0, n_rows * (ROW_SEP / 2.0) + 3.5)

    rearr_rows = classify_chain_rearrangements(chain_rows)
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    ax.set_xlim(-0.32, n_cells_max * (CELL_W + CELL_GAP) + 0.2)
    ax.set_ylim(DB_LABEL_Y_OFFSET - 1.0, n_rows * ROW_SEP + 2.5)
    ax.axis('off')

    db_summary = (f"{len(db_pairs)} deletion bridge pairs (see companion _final_rearr_str_details image)"
                  if len(db_pairs) > 3 else
                  ', '.join(f"DB{d['db_index']}(BP{d['low_bp']}↔BP{d['high_bp']})" for d in db_pairs))
    ax.set_title(
        f'{patient_id} — Baca Chain {chain_number} — FINAL STRUCTURE (after all rearrangements)\n'
        f'n_breakpoints={len(chain_rows)}   n_rearrangements={len(rearr_rows)}   {db_summary}\n'
        f"GAP = native material relocated elsewhere.   Apostrophe = flipped vs. Baca's chain.   "
        f"◆ = real fusion junction.   Purple dotted arc = deletion bridge, current positions.",
        fontsize=13, fontweight='bold', pad=18, loc='center',
    )

    node_positions = {}
    for i, chrom in enumerate(sorted_chroms):
        y = (n_rows - 1 - i) * ROW_SEP
        ax.text(-0.24, y, f'chr{chrom}', ha='right', va='center', fontsize=13, fontweight='bold')
        draw_chromosome_row(ax, y, chrom, cells_by_chrom[chrom], db_lookup, x_start=0.0, node_positions=node_positions)

    draw_deletion_bridge_arcs(ax, db_pairs, node_positions)

    fig.legend(handles=[
        mpatches.Patch(color='#BDBDBD', label='Intact reference segment (unmoved or in-place flip)'),
        mpatches.Patch(color=TELOMERE_COLOR, label='Telomere-ward segment (real, unmapped beyond this point)'),
        mpatches.Patch(facecolor=GAP_COLOR, edgecolor=GAP_HATCH_COLOR, hatch='////', label='GAP — native material relocated'),
        mpatches.Patch(facecolor='white', edgecolor=COL_JUNCTION, label='◆ = real fusion junction'),
        mpatches.Patch(color=DB_PALETTE[0], label='DB-span segment (color = DB pair, same as Part 1)'),
        mpatches.Patch(facecolor='white', edgecolor=DELETION_BRIDGE_COLOR, label='- - - deletion bridge (current positions)'),
    ], loc='lower center', ncol=3, fontsize=10.5, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.05, right=0.98, top=0.83, bottom=0.10)
    out_path = os.path.join(save_dir, f'{patient_id}_chain_{chain_number}_final_rearranged_structure.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


def _add_rearr_type_table(fig, rearr_rows, left, width, top, height, fontsize):
    ax_table = fig.add_axes([left, top - height, width, height])
    ax_table.axis('off')
    ax_table.set_title(f'Rearrangement Type  (n={len(rearr_rows)})', fontsize=12, fontweight='bold', loc='left', pad=6)
    lines = [f"#{r['rearr_number']:<5} BP{r['bp_a']}-BP{r['bp_b']:<6} {r['our_type']}" for r in rearr_rows]
    ax_table.text(0, 1.0, '\n'.join(lines), transform=ax_table.transAxes, ha='left', va='top',
                  fontsize=fontsize, family='monospace', linespacing=1.9)


def _add_segment_table(fig, rows, left, width, top, height, fontsize):
    ax_table = fig.add_axes([left, top - height, width, height])
    ax_table.axis('off')
    ax_table.set_title(f'Segment Spans (native)  (n={len(rows)})', fontsize=12, fontweight='bold', loc='left', pad=6)
    lines = [f'{label:<10} {span}' for label, span in rows]
    ax_table.text(0, 1.0, '\n'.join(lines), transform=ax_table.transAxes, ha='left', va='top',
                  fontsize=fontsize, family='monospace', linespacing=1.9)


def _add_bp_table_wide(fig, chain_rows, left, width, top, height, fontsize):
    rows = chain_rows[['Breakpoint number', 'chromosome', 'position', 'node_strand']].copy()
    rows['Breakpoint number'] = rows['Breakpoint number'].astype(int)
    rows = rows.sort_values('Breakpoint number')
    lines = [f"BP{int(r['Breakpoint number']):<5} chr{int(r['chromosome'])}:{r['position']:,.0f} ({r['node_strand']})"
             for _, r in rows.iterrows()]
    ax_table = fig.add_axes([left, top - height, width, height])
    ax_table.axis('off')
    ax_table.set_title(f'BP # → Location  (n={len(lines)})', fontsize=12, fontweight='bold', loc='left', pad=6)
    ax_table.text(0, 1.0, '\n'.join(lines), transform=ax_table.transAxes, ha='left', va='top',
                  fontsize=fontsize, family='monospace', linespacing=1.9)


def draw_structure_details(patient_id, chain_number, chain_rows, db_pairs, save_dir):
    """Image 2 of the matched pair: the 3 corner tables alone (BP
    location, Segment Spans, Rearrangement Type), laid out as full
    readable columns instead of squeezed corner boxes."""
    rearr_rows = classify_chain_rearrangements(chain_rows)
    use_prefix = len(db_pairs) > 1
    seg_rows = [(f"DB{d['db_index']}:{s['letter']}" if use_prefix else s['letter'],
                 f"spans {s['low_bp']}[+]----{s['high_bp']}[-]")
                for d in db_pairs for s in d['segments']]

    n_bp = len(chain_rows)
    fig_h = max(6.0, max(n_bp, len(rearr_rows), len(seg_rows)) * 0.28 + 2.5)
    fig_w = 18.0

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
    fig.suptitle(f'{patient_id} — Baca Chain {chain_number} — FINAL STRUCTURE DETAILS\n'
                 f'(companion tables for the _final_rearranged_structure image)',
                 fontsize=13, fontweight='bold')

    top, height = 0.88, 0.84
    _add_bp_table_wide(fig, chain_rows, left=0.04, width=0.28, top=top, height=height, fontsize=10.5)
    _add_segment_table(fig, seg_rows, left=0.37, width=0.28, top=top, height=height, fontsize=10.5)
    _add_rearr_type_table(fig, rearr_rows, left=0.70, width=0.28, top=top, height=height, fontsize=10.5)

    out_path = os.path.join(save_dir, f'{patient_id}_chain_{chain_number}_final_rearr_str_details.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


def draw_final_structure_pair(patient_id, chain_number, chain_rows, save_dir):
    db_pairs = get_all_db_pairs_with_segments(chain_rows)
    if not db_pairs:
        return None, None
    db_lookup = _build_db_label_lookup(db_pairs)
    result = assemble_final_structure(chain_rows)
    p1 = draw_structure_diagram(patient_id, chain_number, chain_rows, db_pairs, db_lookup, result, save_dir)
    p2 = draw_structure_details(patient_id, chain_number, chain_rows, db_pairs, save_dir)
    return p1, p2


def run_all(save_dir=OUTPUT_DIR):
    os.makedirs(save_dir, exist_ok=True)
    df = load_mmc5()

    n_drawn, n_skipped_no_db, n_failed = 0, 0, 0
    for (patient_id, chain_number), grp in df.groupby(['Individual', 'Chain number']):
        if not chain_has_db(grp):
            n_skipped_no_db += 1
            continue
        try:
            p1, p2 = draw_final_structure_pair(patient_id, chain_number, grp, save_dir)
            n_drawn += 1
            print(f'{patient_id} chain {chain_number} -> {p1} + {p2}')
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    print(f'\nDrawn: {n_drawn} pairs  Skipped (no DB): {n_skipped_no_db}  Failed: {n_failed}  -> {save_dir}')


if __name__ == '__main__':
    run_all()

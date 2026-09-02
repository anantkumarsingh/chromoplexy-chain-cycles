"""
baca_aberration_cycle_drawing.py  —  AMG cycle structure visualization

Produces TWO output images per patient:

  Image 1  (_observed.png):
    Panel 1 – Full ERG chain (all arcs, all chromosomes)
    Panel 2 – Observed cycle path
    Panel 3 – Observed open path

  Image 2  (_structures.png):
    Panel 1 – Observed open path (detailed)
    Panels 2+ – One schematic per unique cycle structure from enumeration

All panels use SCHEMATIC (equal) spacing so that DSB gaps and arcs are
always visible, even when genomic positions cluster tightly.
"""

import sys, os, math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from baca_aberration_clinical_analysis import (
    load_data,
    _get_erg_chain_rows,
    _build_amg_nodes_and_edges,
    _find_amg_cycles,
    _all_perfect_matchings,
    _cycle_struct_from_matching,
    BACA_CHROM_ABER_FILE_PATH,
    BACA_CLINICAL_PHENO_FILE_PATH,
)

# ── Constants ────────────────────────────────────────────────────────────────

RESULTS_DIR = (
    '/Users/anantkumarsingh/projects/prostate_cancer'
    '/nih-tcga-prad/results/cycle_diagrams'
)

COL = {
    'chr':         '#BDBDBD',
    'chr_faded':   '#E0E0E0',
    'bp':          '#C62828',
    'bp_faded':    '#C4C4C4',
    'telomere':    '#757575',
    'rearr_intra': '#1565C0',
    'rearr_inter': '#E65100',
    'cycle_rearr': '#F57F17',
    'cycle_ref':   '#2E7D32',
    'open_rearr':  '#6A1B9A',
    'open_ref':    '#00838F',
    'open_tel':    '#B71C1C',
    'hyp_cycle':   '#E65100',
    'hyp_open':    '#AAAAAA',
}

X_LO, X_HI = 0.12, 0.88   # chromosome display x-range (leave room for chr label)
GAP         = 0.016        # half-gap either side of a breakpoint mark
CHR_SEP     = 3.5          # vertical distance between chromosome tracks
MAX_REJOINS = 6            # enumeration cutoff


# ── Coordinate helpers ───────────────────────────────────────────────────────

def build_norm_map(sorted_pos):
    """Equal (schematic) spacing — positions evenly spread [X_LO, X_HI]."""
    n = len(sorted_pos)
    if n == 1:
        return {sorted_pos[0]: (X_LO + X_HI) / 2}
    return {p: X_LO + i * (X_HI - X_LO) / (n - 1)
            for i, p in enumerate(sorted_pos)}


def node_x(pos, strand, nmap):
    return nmap[pos] + (GAP if strand == '+' else -GAP)


# ── Drawing primitives ───────────────────────────────────────────────────────

def _arc(ax, xa, ya, xb, yb, rad, color, lw, ls='--', alpha=0.85,
         arrowstyle='-', zorder=4):
    ax.add_patch(FancyArrowPatch(
        (xa, ya), (xb, yb),
        connectionstyle=f'arc3,rad={rad}',
        arrowstyle=arrowstyle, mutation_scale=15,
        color=color, linewidth=lw, linestyle=ls, alpha=alpha, zorder=zorder,
    ))


def _badge(ax, x, y, text, color):
    """Step-label badge with white rounded box."""
    ax.text(x, y, text, ha='center', va='bottom', fontsize=8, fontweight='bold',
            color=color,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=color, lw=1.1, alpha=0.93))


def _ref_highlight(ax, x_lo, x_hi, y, going_right, color, label):
    """Highlight a chromosome reference segment + direction arrow + badge."""
    ax.plot([x_lo, x_hi], [y, y], color=color, linewidth=11,
            solid_capstyle='butt', zorder=5, alpha=0.65)
    mid, d = (x_lo + x_hi) / 2, 0.035
    ax.annotate('', xy=(mid + (d if going_right else -d), y),
                xytext=(mid + (-d if going_right else d), y),
                arrowprops=dict(arrowstyle='-|>', color=color,
                                lw=2.0, mutation_scale=14), zorder=8)
    _badge(ax, mid, y - 0.52, label, color)


# ── Chromosome track ─────────────────────────────────────────────────────────

def draw_chr_track(ax, chrom, y, sorted_pos, nmap, faded=False):
    bar = COL['chr_faded'] if faded else COL['chr']
    bounds = [X_LO] + [nmap[p] for p in sorted_pos] + [X_HI]
    for i, (a, b) in enumerate(zip(bounds[:-1], bounds[1:])):
        s = a + (GAP if i > 0 else 0)
        e = b - (GAP if i < len(bounds) - 2 else 0)
        ax.plot([s, e], [y, y], color=bar, linewidth=11,
                solid_capstyle='butt', zorder=1)
    for p in sorted_pos:
        ax.plot([nmap[p], nmap[p]], [y - 0.22, y + 0.22],
                color=COL['bp_faded'] if faded else COL['bp'],
                linewidth=2.5, zorder=3)
    for xe in (X_LO, X_HI):
        ax.plot([xe, xe], [y - 0.15, y + 0.15],
                color=COL['telomere'], linewidth=4, solid_capstyle='round', zorder=2)
    ax.text(X_LO - 0.07, y, f'chr{int(chrom)}', ha='right', va='center',
            fontsize=10, fontweight='bold',
            color='#AAAAAA' if faded else '#222222')


def draw_bp_labels(ax, chrom, y, sorted_pos, nmap, rearr_adj, faded=False, max_positions=40):
    """
    Staggered labels: even-indexed DSBs labelled ABOVE, odd-indexed BELOW.
    Skips entirely past max_positions — dense chromosomes (genome-wide
    patients can have 100+ DSBs on one chromosome) would render an
    unreadable, very slow label pile; the tick marks from draw_chr_track
    still show DSB density.
    """
    if len(sorted_pos) > max_positions:
        return
    txt = '#C0C0C0' if faded else '#444444'
    bpc = '#BBBBBB' if faded else COL['bp']
    for i, pos in enumerate(sorted_pos):
        x    = nmap[pos]
        above = (i % 2 == 0)
        ly   = y + 0.38 if above else y - 0.38
        va   = 'bottom'  if above else 'top'
        ax.text(x, ly, f'{pos/1e6:.2f} Mb', ha='center', va=va,
                fontsize=6.5, color=txt, rotation=0)
        strands = [s for s in ('+', '-') if (chrom, pos, s) in rearr_adj]
        if strands and not faded:
            sy = y + 0.18 if above else y - 0.18
            sva = 'bottom' if above else 'top'
            ax.text(x, sy, '(' + '/'.join(strands) + ')',
                    ha='center', va=sva, fontsize=7, color=bpc, fontweight='bold')


def _faded_backdrop(ax, sorted_chroms, chr_y, chr_to_dsb_pos, nmaps, rearr_adj):
    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c], faded=True)
        draw_bp_labels(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c], rearr_adj, faded=True)


# ── Arc curvature helpers ─────────────────────────────────────────────────────

def _assign_rads(arc_list, is_inter_list,
                 rad_min=-0.22, rad_max=-0.78):
    """
    Sort arcs by x-span, assign curvature shallow→deep.
    Returns list of rad values parallel to arc_list.
    """
    n = len(arc_list)
    intra_indices = [i for i, inter in enumerate(is_inter_list) if not inter]
    # Sort intra arcs by span
    intra_sorted = sorted(intra_indices,
                          key=lambda i: abs(arc_list[i][2] - arc_list[i][0]))
    rads = [0.28] * n   # default for inter
    for rank, idx in enumerate(intra_sorted):
        rads[idx] = (rad_min + rank * (rad_max - rad_min) / max(len(intra_sorted) - 1, 1))
    return rads


# ── Observed cycle path drawing ───────────────────────────────────────────────

def _draw_cycle_hop(ax, n_from, n_to, hop_type, step_n, nmaps, chr_y):
    if hop_type == 'rearr':
        xa = node_x(n_from[1], n_from[2], nmaps[n_from[0]])
        xb = node_x(n_to[1],   n_to[2],   nmaps[n_to[0]])
        ya, yb = chr_y[n_from[0]], chr_y[n_to[0]]
        inter = abs(ya - yb) > 0.1
        rad = 0.28 if inter else -0.55
        _arc(ax, xa, ya, xb, yb, rad, COL['cycle_rearr'], 2.6, '-', 1.0, '-|>', 7)
        _badge(ax, (xa + xb) / 2,
               (ya + yb) / 2 + (0.62 if not inter else 0.30),
               f'[{step_n}] REARR', COL['cycle_rearr'])
    else:
        lo, hi = min(n_from[1], n_to[1]), max(n_from[1], n_to[1])
        nmap   = nmaps[n_from[0]]
        y      = chr_y[n_from[0]]
        _ref_highlight(ax, nmap[lo] + GAP, nmap[hi] - GAP, y,
                       n_to[1] > n_from[1], COL['cycle_ref'], f'[{step_n}] REF')


def _draw_open_hop(ax, n_from, n_to, hop_type, step_n, nmaps, chr_y):
    if hop_type == 'rearr':
        xa = node_x(n_from[1], n_from[2], nmaps[n_from[0]])
        xb = node_x(n_to[1],   n_to[2],   nmaps[n_to[0]])
        ya, yb = chr_y[n_from[0]], chr_y[n_to[0]]
        inter  = abs(ya - yb) > 0.1
        rad    = 0.28 if inter else -0.55
        _arc(ax, xa, ya, xb, yb, rad, COL['open_rearr'], 2.6, '-', 1.0, '-|>', 7)
        _badge(ax, (xa + xb) / 2,
               (ya + yb) / 2 + (0.62 if not inter else 0.30),
               f'[{step_n}] REARR', COL['open_rearr'])
    else:
        lo, hi = min(n_from[1], n_to[1]), max(n_from[1], n_to[1])
        nmap   = nmaps[n_from[0]]
        y      = chr_y[n_from[0]]
        _ref_highlight(ax, nmap[lo] + GAP, nmap[hi] - GAP, y,
                       n_to[1] > n_from[1], COL['open_ref'], f'[{step_n}] REF')


def _draw_telomere(ax, last_node, side, nmaps, chr_y, step_n):
    chrom, pos, strand = last_node
    x0    = node_x(pos, strand, nmaps[chrom])
    y     = chr_y[chrom]
    x_end = X_HI + 0.08 if side == 'right' else X_LO - 0.08
    mid   = (x0 + x_end) / 2
    _arc(ax, x0, y, x_end, y, 0.0, COL['open_tel'], 1.8, '--', 0.9, '-|>', 7)
    ax.text(x_end, y, '⊣' if side == 'right' else '⊢',
            ha='left' if side == 'right' else 'right',
            va='center', fontsize=16, color=COL['open_tel'], fontweight='bold', zorder=9)
    _badge(ax, mid, y + 0.44, f'[{step_n}] REF→TEL', COL['open_tel'])


def _start_marker(ax, node, nmaps, chr_y, color):
    sx = node_x(node[1], node[2], nmaps[node[0]])
    sy = chr_y[node[0]]
    ax.plot(sx, sy, 'o', color=color, markersize=11,
            markeredgecolor='white', markeredgewidth=2, zorder=9)
    ax.text(sx, sy + 0.50, 'START', ha='center', va='bottom',
            fontsize=7.5, fontweight='bold', color=color)


# ── Unified traversal ─────────────────────────────────────────────────────────

def extract_all_traversal_hops(all_nodes, rearr_adj, ref_adj):
    visited, cycle_hops, open_paths = set(), [], []
    for start in sorted(all_nodes):
        if start in visited or start not in rearr_adj:
            continue
        path, hops = [], []
        current, htype, rearr_hops = start, 'rearrangement', 0
        is_open = termination = telomere_side = None
        while current not in visited:
            visited.add(current); path.append(current)
            if htype == 'rearrangement':
                nxt = rearr_adj.get(current)
                if nxt is None:
                    is_open = True; termination = 'MISSING_REARR'; break
                hops.append((current, nxt, 'rearr'))
                rearr_hops += 1; htype = 'reference'
            else:
                nxt = ref_adj.get(current)
                if nxt is None or nxt == 'TELOMERE':
                    is_open = True; termination = 'TELOMERE'
                    telomere_side = 'right' if current[2] == '+' else 'left'; break
                hops.append((current, nxt, 'ref')); htype = 'rearrangement'
            current = nxt
        if not is_open and current == start and rearr_hops > 0:
            cycle_hops.append(hops)
        elif is_open:
            open_paths.append({'hops': hops, 'last_node': current,
                                'termination': termination,
                                'telomere_side': telomere_side})
    return cycle_hops, open_paths


# ── Enumeration ───────────────────────────────────────────────────────────────

def find_one_example_per_structure(rearr_adj, all_nodes, ref_adj):
    endpoints = list(rearr_adj.keys())
    n_rejoins = len(endpoints) // 2
    if n_rejoins > MAX_REJOINS:
        return None, None, None
    examples, counts, total = {}, defaultdict(int), 0
    for matching in _all_perfect_matchings(endpoints):
        struct = _cycle_struct_from_matching(matching, ref_adj, all_nodes)
        counts[struct] += 1; total += 1
        if struct not in examples:
            examples[struct] = matching
    return examples, counts, total


def cycle_arcs_in_matching(matching, all_nodes, ref_adj):
    temp = {u: v for u, v in matching}
    temp.update({v: u for u, v in matching})
    cycle_hops, _ = extract_all_traversal_hops(all_nodes, temp, ref_adj)
    pairs = set()
    for hops in cycle_hops:
        for nf, nt, ht in hops:
            if ht == 'rearr':
                pairs.add((nf, nt)); pairs.add((nt, nf))
    return pairs


# ── Panel drawing ─────────────────────────────────────────────────────────────

def _ax_defaults(ax, y_span):
    ax.set_facecolor('white')
    ax.set_xlim(-0.02, 1.10)
    ax.set_ylim(-1.8 - y_span * 0.05, y_span + 2.2)
    ax.axis('off')


def draw_full_chain(ax, sorted_chroms, chr_y, chr_to_dsb_pos, nmaps,
                    rearr_adj, cycle_str, n_dsbs, rank, n_open, theta,
                    title_prefix='Full ERG Chain'):
    ax.set_title(f'{title_prefix}\n{theta}\nn={n_dsbs}  Rank={rank}'
                 f'  Observed: {cycle_str}  ({n_open} open)',
                 fontsize=9, fontweight='bold', pad=12)
    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c])
        draw_bp_labels(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c], rearr_adj)

    # Collect arcs and assign curvature by span
    seen, arc_data = set(), []
    for u, v in rearr_adj.items():
        key = tuple(sorted([u, v]))
        if key in seen: continue
        seen.add(key)
        xa = node_x(u[1], u[2], nmaps[u[0]])
        xb = node_x(v[1], v[2], nmaps[v[0]])
        ya, yb = chr_y[u[0]], chr_y[v[0]]
        arc_data.append((xa, ya, xb, yb, u[0] != v[0], u, v))

    is_inter_list = [a[4] for a in arc_data]
    rads = _assign_rads(arc_data, is_inter_list)

    for (xa, ya, xb, yb, is_inter, u, v), rad in zip(arc_data, rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 1.8, '--', 0.85, '-', 4)

    ax.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Interchromosomal'),
    ], loc='lower right', fontsize=8, framealpha=0.95, edgecolor='#CCCCCC')


def draw_observed_cycle(ax, cycle_hops, sorted_chroms, chr_y, chr_to_dsb_pos,
                        nmaps, rearr_adj, cycle_str):
    ax.set_title(f'Observed — Cycle Path\n{cycle_str}\n'
                 f'Amber = REARR  ·  Green = REF',
                 fontsize=9, fontweight='bold', pad=12)
    _faded_backdrop(ax, sorted_chroms, chr_y, chr_to_dsb_pos, nmaps, rearr_adj)
    if not cycle_hops:
        ax.text(0.5, 0.5, 'No closed cycles\nin observed data',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=11, color='#888888', style='italic')
        return
    for hops in cycle_hops:
        for sn, (nf, nt, ht) in enumerate(hops, 1):
            _draw_cycle_hop(ax, nf, nt, ht, sn, nmaps, chr_y)
        _start_marker(ax, hops[0][0], nmaps, chr_y, COL['cycle_rearr'])
    ax.legend(handles=[
        mpatches.Patch(color=COL['cycle_rearr'], label='REARR hop'),
        mpatches.Patch(color=COL['cycle_ref'],   label='REF segment'),
    ], loc='lower right', fontsize=8, framealpha=0.95, edgecolor='#CCCCCC')


def draw_observed_open(ax, open_path_data, sorted_chroms, chr_y,
                       chr_to_dsb_pos, nmaps, rearr_adj, n_open):
    ax.set_title(f'Observed — Open Path  (1 of {n_open})\n'
                 f'Purple = REARR  ·  Teal = REF  ·  Red ⊣ = TEL',
                 fontsize=9, fontweight='bold', pad=12)
    _faded_backdrop(ax, sorted_chroms, chr_y, chr_to_dsb_pos, nmaps, rearr_adj)
    if not open_path_data:
        ax.text(0.5, 0.5, 'No open paths', ha='center', va='center',
                transform=ax.transAxes, fontsize=11, color='#888888', style='italic')
        return
    op   = open_path_data[0]
    hops = op['hops']
    for sn, (nf, nt, ht) in enumerate(hops, 1):
        _draw_open_hop(ax, nf, nt, ht, sn, nmaps, chr_y)
    next_step = len(hops) + 1
    if op['termination'] == 'TELOMERE':
        _draw_telomere(ax, op['last_node'], op['telomere_side'],
                       nmaps, chr_y, next_step)
    if hops:
        _start_marker(ax, hops[0][0], nmaps, chr_y, COL['open_rearr'])
    ax.legend(handles=[
        mpatches.Patch(color=COL['open_rearr'], label='REARR hop'),
        mpatches.Patch(color=COL['open_ref'],   label='REF segment'),
        mpatches.Patch(color=COL['open_tel'],   label='TELOMERE'),
    ], loc='lower right', fontsize=8, framealpha=0.95, edgecolor='#CCCCCC')


def draw_structure_panel(ax, matching, struct_tuple, count, total,
                         sorted_chroms, chr_y, chr_to_dsb_pos,
                         nmaps, rearr_adj, all_nodes, ref_adj):
    label = (' + '.join(f'C{n}' for n in struct_tuple)
             if struct_tuple else 'none')
    pct = count / total * 100
    ax.set_title(f'{label}\n{count}/{total}  ({pct:.1f}%)',
                 fontsize=9, fontweight='bold', pad=10)
    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_dsb_pos[c], nmaps[c], faded=True)

    cyc_pairs = cycle_arcs_in_matching(matching, all_nodes, ref_adj)

    # Build arc list with display coords
    arc_data, is_inter_list = [], []
    for u, v in matching:
        xa = node_x(u[1], u[2], nmaps[u[0]])
        xb = node_x(v[1], v[2], nmaps[v[0]])
        ya, yb = chr_y[u[0]], chr_y[v[0]]
        inter = u[0] != v[0]
        arc_data.append((xa, ya, xb, yb, inter, u, v))
        is_inter_list.append(inter)

    rads = _assign_rads(arc_data, is_inter_list)

    for (xa, ya, xb, yb, inter, u, v), rad in zip(arc_data, rads):
        is_cycle = (u, v) in cyc_pairs
        col  = COL['hyp_cycle'] if is_cycle else COL['hyp_open']
        lw   = 2.4              if is_cycle else 1.4
        ls   = '-'              if is_cycle else '--'
        al   = 0.88             if is_cycle else 0.50
        _arc(ax, xa, ya, xb, yb, rad, col, lw, ls, al, '-', 5)

    ax.legend(handles=[
        mpatches.Patch(color=COL['hyp_cycle'], label='cycle arc'),
        mpatches.Patch(color=COL['hyp_open'],  label='open arc'),
    ], loc='lower right', fontsize=7, framealpha=0.9,
       edgecolor='#CCCCCC', handlelength=1.2)


# ── Main drawing ──────────────────────────────────────────────────────────────

def draw_patient_amg(patient_id, chrom_df, clinical_df, save_dir=RESULTS_DIR):
    os.makedirs(save_dir, exist_ok=True)

    # ── AMG computation ───────────────────────────────────────────────────────
    chain_rows = _get_erg_chain_rows(patient_id, chrom_df, clinical_df)
    if chain_rows.empty:
        print(f'[{patient_id}] No chain rows.'); return

    all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(chain_rows)
    cycles, _         = _find_amg_cycles(all_nodes, rearr_adj, ref_adj)
    cycle_hops, opaths = extract_all_traversal_hops(all_nodes, rearr_adj, ref_adj)
    examples, counts, total = find_one_example_per_structure(
        rearr_adj, all_nodes, ref_adj)

    n_rejoins = len(chain_rows)
    n_dsbs    = n_rejoins * 2
    k_local   = len(chr_to_dsb_pos)
    b_i       = {c: len(v) for c, v in chr_to_dsb_pos.items()}
    cycle_str = (' + '.join(f'C{n}' for n in sorted(cycles, reverse=True))
                 or 'none')
    theta     = (f"Θ({k_local}, "
                 f"({', '.join(str(b_i[c]) for c in sorted(chr_to_dsb_pos))}))")
    rank      = str(n_dsbs - len(cycles)) if cycles else 'N/A'

    # ── Layout helpers ────────────────────────────────────────────────────────
    sorted_chroms = sorted(chr_to_dsb_pos)
    chr_y  = {c: (len(sorted_chroms) - 1 - i) * CHR_SEP
              for i, c in enumerate(sorted_chroms)}
    nmaps  = {c: build_norm_map(chr_to_dsb_pos[c]) for c in sorted_chroms}
    y_span = (len(sorted_chroms) - 1) * CHR_SEP

    panel_w = 8.0          # width per panel (generous for labels)
    panel_h = max(5.0, k_local * 3.2 + 2.5)   # height scales with #chromosomes

    # ── IMAGE 1: 3 fixed panels ───────────────────────────────────────────────
    fig1, axes1 = plt.subplots(1, 3, figsize=(3 * panel_w, panel_h))
    for ax in axes1:
        _ax_defaults(ax, y_span)

    draw_full_chain(axes1[0], sorted_chroms, chr_y, chr_to_dsb_pos, nmaps,
                    rearr_adj, cycle_str, n_dsbs, rank, len(opaths), theta)
    draw_observed_cycle(axes1[1], cycle_hops, sorted_chroms, chr_y,
                        chr_to_dsb_pos, nmaps, rearr_adj, cycle_str)
    draw_observed_open(axes1[2], opaths, sorted_chroms, chr_y,
                       chr_to_dsb_pos, nmaps, rearr_adj, len(opaths))

    fig1.suptitle(
        f'{patient_id}  —  Observed AMG  |  {theta}  |  n={n_dsbs}',
        fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout(pad=2.0)
    out1 = os.path.join(save_dir, f'{patient_id}_observed.png')
    fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig1)
    print(f'Image 1 → {out1}')

    # ── IMAGE 2: open path + one panel per structure ──────────────────────────
    if examples is None:
        print(f'[{patient_id}] Skipping Image 2 (n_rejoins={n_rejoins} > {MAX_REJOINS}).'); return

    def _struct_key(s):
        return (sum(s), -counts.get(s, 0))

    sorted_structs = sorted(examples.keys(), key=_struct_key)
    n_structs = len(sorted_structs)
    n_panels  = 1 + n_structs      # open path + each structure
    n_cols    = min(4, n_panels)
    n_rows    = math.ceil(n_panels / n_cols)

    struct_panel_w = 5.5
    struct_panel_h = max(4.5, k_local * 2.8 + 2.0)

    fig2, axes_grid = plt.subplots(n_rows, n_cols,
                                   figsize=(n_cols * struct_panel_w,
                                            n_rows * struct_panel_h),
                                   squeeze=False)
    axes2 = [axes_grid[r][c] for r in range(n_rows) for c in range(n_cols)]
    for ax in axes2:
        _ax_defaults(ax, y_span)

    # Panel 0: observed open path
    draw_observed_open(axes2[0], opaths, sorted_chroms, chr_y,
                       chr_to_dsb_pos, nmaps, rearr_adj, len(opaths))

    # Panels 1+: one per unique cycle structure
    for i, struct in enumerate(sorted_structs):
        draw_structure_panel(
            axes2[1 + i], examples[struct], struct, counts[struct], total,
            sorted_chroms, chr_y, chr_to_dsb_pos,
            nmaps, rearr_adj, all_nodes, ref_adj,
        )

    # Hide unused axes
    for ax in axes2[n_panels:]:
        ax.set_visible(False)

    fig2.suptitle(
        f'{patient_id}  —  Cycle Structure Atlas  |  {theta}  |  '
        f'Total matchings: {total}',
        fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout(pad=2.0)
    out2 = os.path.join(save_dir, f'{patient_id}_structures.png')
    fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig2)
    print(f'Image 2 → {out2}  ({n_structs} unique structures)')


# ── Batch runner ─────────────────────────────────────────────────────────────

def run_all_patients(chrom_df, clinical_df, save_dir=RESULTS_DIR):
    """
    Run draw_patient_amg for every patient that has a non-empty ERG chain.
    Patients with n_rejoins > MAX_REJOINS produce Image 1 only (no enumeration).
    """
    all_ids = sorted(chrom_df['Individual'].unique(), key=str)
    # Pre-filter to those with a non-empty chain so we know total upfront
    has_chain = []
    for pid in all_ids:
        chain = _get_erg_chain_rows(pid, chrom_df, clinical_df)
        if not chain.empty:
            has_chain.append((pid, len(chain)))

    n_total    = len(has_chain)
    succeeded  = []
    img2_done  = []
    img2_skip  = []
    failed     = []

    print(f"\n{'='*60}")
    print(f"  AMG Cycle Drawing  —  {n_total} patients with ERG chains")
    print(f"{'='*60}\n")

    for i, (pid, nr) in enumerate(has_chain, 1):
        # Double factorial for display
        m = 1
        for k in range(1, 2 * nr, 2):
            m *= k
        enumerable = nr <= MAX_REJOINS
        tag = f"{m:,} matchings" if enumerable else f"n_rejoins={nr} (Image 2 skipped)"
        print(f"[{i:>2}/{n_total}] {pid:<30}  {tag}")

        try:
            draw_patient_amg(pid, chrom_df, clinical_df, save_dir=save_dir)
            succeeded.append(pid)
            if enumerable:
                img2_done.append(pid)
            else:
                img2_skip.append(pid)
        except Exception as exc:
            print(f"         ERROR: {exc}")
            failed.append((pid, str(exc)))

    print(f"\n{'='*60}")
    print(f"  Done.")
    print(f"  Succeeded : {len(succeeded)}/{n_total}")
    print(f"  Image 2   : {len(img2_done)} patients (enumerable)")
    print(f"  Image 1   : {len(img2_skip)} patients (too large to enumerate)")
    if failed:
        print(f"  Failed    : {len(failed)}")
        for pid, err in failed:
            print(f"    {pid}: {err}")
    print(f"{'='*60}\n")
    return succeeded, failed


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    chrom_df, clinical_df = load_data(
        BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)
    run_all_patients(chrom_df, clinical_df)

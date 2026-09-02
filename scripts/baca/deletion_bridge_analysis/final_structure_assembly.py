"""
final_structure_assembly.py — Part 2's core algorithm: build the final
post-rearrangement derivative structure for a whole mmc5 chain, by
applying that chain's real rearrangement edges in ASCENDING
`Rearrangement number` order (locked build-order rule, confirmed
2026-08-13 — see DELETION_BRIDGE_ANALYSIS_PLAN.md decision #4).

Why this is a NEW algorithm, not a reuse of
core/baca_aberration_cycle_drawing.py's extract_all_traversal_hops:
that function treats every node with a real rearrangement partner as an
INDEPENDENT traversal starting point, and never traverses a node's own
reference edge once it has already started a trace from it. This means a
segment whose BOTH ends independently carry real rearrangement partners
(e.g. Piece 1 in P05-1657: 2+ fuses to 14-, 3- fuses to 15-) gets
reported as belonging to TWO SEPARATE disconnected fragments instead of
being recognized as the single physical bridge connecting them —
verified directly by running that function on P05-1657 this session
(CONTEXT_LOG_003.md Session 8). This module instead ASSEMBLES molecules
by literally simulating breakage-and-rejoining: every reference segment
starts as its own molecule, and each real rearrangement edge (processed
in the fixed order above) merges two molecules at their matching free
ends, correctly re-orienting (flipping) one side when needed. Because
every segment's own body is a first-class object being merged (not just
node-to-node graph traversal), a "both ends real-fused" segment is
naturally captured as an internal, correctly-oriented piece of the
merged molecule — this is what fixes the Piece 1 problem.

Nodes are represented as (breakpoint_number, strand) tuples -- simpler
than (chrom, pos, strand) since mmc5 breakpoint numbers are already
unique within a chain and already encode chromosome+position via the
chain's own breakpoint table.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from db_segments_common import load_mmc5  # noqa: E402


class Piece:
    """One intact reference segment (native = low_bp+ ---- high_bp-),
    possibly flipped (high_bp- ---- low_bp+) after assembly.

    `moved` (added 2026-08-13, per the user's exact generative process for
    the P05-1657 example): True if, at ANY point during assembly, this
    piece's unit was the one that had to relocate to reach its final
    connection (as opposed to a "local" merge with an already-native
    neighbor, which never sets this). Drives the GAP-vs-shown-in-place
    decision for Part 2's chromosome-row visualization. Set once, never
    cleared, and carried through flip()."""
    __slots__ = ('chrom', 'low_bp', 'high_bp', 'flipped', 'moved')

    def __init__(self, chrom, low_bp, high_bp, flipped=False, moved=False):
        self.chrom = chrom
        self.low_bp = low_bp
        self.high_bp = high_bp
        self.flipped = flipped
        self.moved = moved

    @property
    def left_node(self):
        return (self.high_bp, '-') if self.flipped else (self.low_bp, '+')

    @property
    def right_node(self):
        return (self.low_bp, '+') if self.flipped else (self.high_bp, '-')

    def display_nodes(self, display_swapped):
        """(left, right) boundary nodes to draw for THIS specific row
        cell. `display_swapped` is a pure rendering-direction choice
        (True when this piece sits in a molecule being entered from its
        CHAIN_TERMINAL side on a different chromosome's row -- see
        build_chromosome_rows) and is completely independent of
        `.flipped` (which stays fixed, driving the apostrophe/DB-lookup
        semantics everywhere else). Swapping display order here never
        mutates the piece or touches `.flipped`."""
        return (self.right_node, self.left_node) if display_swapped else (self.left_node, self.right_node)

    def flip(self):
        return Piece(self.chrom, self.low_bp, self.high_bp, not self.flipped, self.moved)

    def native_span_text(self):
        return f'{self.low_bp}[+]----{self.high_bp}[-]'

    def current_span_text(self):
        if self.flipped:
            return f'{self.high_bp}[-]----{self.low_bp}[+]'
        return f'{self.low_bp}[+]----{self.high_bp}[-]'

    def __repr__(self):
        return f'Piece(chr{self.chrom}, {self.low_bp}-{self.high_bp}, flipped={self.flipped})'


class Molecule:
    """Ordered left-to-right list of Pieces (possibly empty -- a bare,
    unconsumed telomere-facing tip), with its two current free ends."""
    __slots__ = ('pieces', 'left_free', 'right_free')

    def __init__(self, pieces, left_free, right_free):
        self.pieces = pieces
        self.left_free = left_free
        self.right_free = right_free

    def reversed(self):
        return Molecule([p.flip() for p in reversed(self.pieces)], self.right_free, self.left_free)


def build_initial_molecules(chain_rows):
    """One molecule per intact reference segment (position-adjacent
    breakpoints, same chromosome, native/unflipped orientation), plus a
    trivial (0-piece) molecule for every genuinely telomere-facing tip
    (a chromosome's lowest breakpoint's '-' side, highest breakpoint's
    '+' side within this chain) so a real fusion edge landing there still
    has a molecule to merge into.

    Returns: dict node -> Molecule (both of a molecule's free ends map to
    the same Molecule object; multi-piece molecules initially all have
    exactly 1 piece).

    IMPORTANT (bug found and fixed 2026-08-13, caught running the full
    194-chain cohort, not just the P05-1657 test case): breakpoint
    POSITION is not a safe dedup key. Two distinct breakpoint numbers can
    share the exact same (chromosome, position) -- a known data quirk
    already documented in CLAUDE.md's "Strand convention verified"
    section (e.g. P08-1042 chr2:65,557,230 is breakpoint 24 for
    Rearrangement 31 AND breakpoint 77 for Rearrangement 53). Grouping by
    `position.unique()` silently drops one of the two breakpoints, which
    then crashes assembly when a later rearrangement needs the dropped
    one. Fixed by sorting each chromosome's breakpoints by
    (position, breakpoint_number) directly -- every distinct breakpoint
    number is preserved even when two share a position (producing a
    harmless zero-length segment between them).
    """
    free_end_to_molecule = {}
    for chrom, grp in chain_rows.groupby('chromosome'):
        chrom = int(chrom)
        bps_sorted = [int(bp) for bp in
                      grp[['Breakpoint number', 'position']]
                      .drop_duplicates(subset=['Breakpoint number'])
                      .sort_values(['position', 'Breakpoint number'])['Breakpoint number']]
        for i in range(len(bps_sorted) - 1):
            low_bp, high_bp = bps_sorted[i], bps_sorted[i + 1]
            piece = Piece(chrom, low_bp, high_bp, flipped=False)
            mol = Molecule([piece], piece.left_node, piece.right_node)
            free_end_to_molecule[piece.left_node] = mol
            free_end_to_molecule[piece.right_node] = mol
        lowest_bp, highest_bp = bps_sorted[0], bps_sorted[-1]
        for tip_node in [(lowest_bp, '-'), (highest_bp, '+')]:
            if tip_node not in free_end_to_molecule:
                free_end_to_molecule[tip_node] = Molecule([], tip_node, tip_node)
    return free_end_to_molecule


def build_bp_position_lookup(chain_rows):
    """bp -> (chrom, position), used only for the relocation tiebreak
    (which of two substantive units is "earlier" and therefore tagged
    the mover) -- see DELETION_BRIDGE_ANALYSIS_PLAN.md's Part 2 gap rule,
    confirmed 2026-08-13 against the user's own P05-1657 derivation."""
    bp_position = {}
    for chrom, grp in chain_rows.groupby('chromosome'):
        chrom = int(chrom)
        for bp, pos in zip(grp['Breakpoint number'].astype(int), grp['position']):
            bp_position[bp] = (chrom, pos)
    return bp_position


def _unit_bps(unit):
    """All breakpoint numbers currently present in a unit (Molecule) --
    from its real pieces, or its own single node if trivial."""
    if unit.pieces:
        bps = set()
        for p in unit.pieces:
            bps.add(p.low_bp)
            bps.add(p.high_bp)
        return bps
    return {unit.left_free[0]}  # trivial molecule: left_free == right_free == its own node


def classify_rearrangement_type(same_chrom, same_strand):
    """Independent 4-way classification from chromosome-match + strand-
    match alone -- NOT from Baca's own `Class` field (may be cross-
    checked against it afterward via classify_rearrangements_with_baca_
    crosscheck, but that comparison must never influence this
    classification). Rule confirmed 2026-08-13, reverse-engineered from
    how the user classified P05-1657's 4 real rearrangements by hand (3/4
    cases verified to match Baca's Class exactly; the diff-chrom+
    opposite-strand case was never observed there and is an inferred
    completion of the pattern -- see DELETION_BRIDGE_ANALYSIS_PLAN.md
    decision #7)."""
    if same_chrom and not same_strand:
        return 'SIMPLE / NON-RECIPROCAL TRANSLOCATION'
    if same_chrom and same_strand:
        return 'SIMPLE INVERSION'
    if not same_chrom and same_strand:
        return 'INVERTED TRANSLOCATION'
    return 'TRANSLOCATION'  # different chromosome, opposite strand -- inferred 4th category


def classify_chain_rearrangements(chain_rows):
    """One row per real fusion edge in this chain (any Rearrangement
    number with exactly 2 breakpoints): our independent type + Baca's own
    Class field alongside it for comparison only."""
    rows = []
    for rearr_num, grp in chain_rows.groupby('Rearrangement number'):
        if len(grp) != 2:
            continue
        r1, r2 = grp.iloc[0], grp.iloc[1]
        same_chrom = int(r1['chromosome']) == int(r2['chromosome'])
        same_strand = r1['node_strand'] == r2['node_strand']
        our_type = classify_rearrangement_type(same_chrom, same_strand)
        rows.append({
            'rearr_number': rearr_num,
            'bp_a': int(r1['Breakpoint number']), 'bp_b': int(r2['Breakpoint number']),
            'chrom_a': int(r1['chromosome']), 'chrom_b': int(r2['chromosome']),
            'strand_a': r1['node_strand'], 'strand_b': r2['node_strand'],
            'our_type': our_type,
        })
    rows.sort(key=lambda r: r['rearr_number'])
    return rows


def build_rearr_edges(chain_rows):
    """List of (rearr_number, node_u, node_v) — one per real fusion edge
    in this chain, sorted ASCENDING by rearr_number (locked build order).
    Rearrangements with a partner outside this chain (unresolved, len!=2)
    are skipped -- that end stays a free/open tip, correctly reflecting
    that this chain's real data doesn't resolve it."""
    edges = []
    n_unresolved = 0
    for rearr_num, grp in chain_rows.groupby('Rearrangement number'):
        bps = grp['Breakpoint number'].astype(int).tolist()
        strands = grp['node_strand'].tolist()
        if len(bps) != 2:
            n_unresolved += len(bps)
            continue
        node_a = (bps[0], strands[0])
        node_b = (bps[1], strands[1])
        edges.append((rearr_num, node_a, node_b))
    edges.sort(key=lambda e: e[0])
    return edges, n_unresolved


def assemble_final_structure(chain_rows):
    """
    Runs the full assembly for one chain (all its breakpoints/rearrangements
    -- NOT scoped to any single DB pair, since real material from anywhere
    in the chain can move into a DB span, as verified for Piece 1 in
    P05-1657).

    Returns: {
      'open_molecules': [Molecule, ...],   # final derivative fragments, still open
      'closed_molecules': [(rearr_num, Molecule), ...],  # self-closed into a cycle
      'rearr_log': [(rearr_num, node_u, node_v, event), ...],  # for verification/debugging
      'n_unresolved_rearr_ends': int,
    }
    """
    free_end_to_molecule = build_initial_molecules(chain_rows)
    rearr_edges, n_unresolved = build_rearr_edges(chain_rows)
    bp_position = build_bp_position_lookup(chain_rows)

    closed_molecules = []
    fully_capped_molecules = []
    rearr_log = []

    for rearr_num, node_u, node_v in rearr_edges:
        M = free_end_to_molecule.get(node_u)
        N = free_end_to_molecule.get(node_v)
        if M is None or N is None:
            raise AssertionError(
                f'node already consumed or missing while processing rearr {rearr_num}: '
                f'{node_u} (found={M is not None}) / {node_v} (found={N is not None})'
            )

        if M is N:
            for end in {M.left_free, M.right_free}:
                free_end_to_molecule.pop(end, None)
            closed_molecules.append((rearr_num, M))
            rearr_log.append((rearr_num, node_u, node_v, 'SELF_CLOSURE'))
            continue

        # GAP classification (added 2026-08-13, confirmed against the
        # user's own P05-1657 derivation): "local" if the breakpoint
        # being fused on ONE side is already a boundary of some piece
        # CURRENTLY in the OTHER unit -- meaning that unit already has
        # real structure sitting immediately next to this exact
        # breakpoint, so this fusion is just an in-place flip/rewire
        # between already-adjacent content, not a relocation. (An
        # earlier version of this check computed the NATIVE reference
        # neighbor of the specific strand-side instead -- wrong, and
        # caught by testing: it misclassified Rearr 17 in P05-1657 as a
        # non-local relocation and wrongly flagged Segment C as moved,
        # when C should stay unmarked since B -- already fused to C's
        # own unit via Rearr 11 -- is C's real native neighbor. Checking
        # "is the breakpoint itself already a member of the other unit"
        # directly captures that, with no extra indirection needed.)
        # Otherwise it's a "relocation" event: the earlier-native-
        # position unit is tagged as the mover (its pieces get
        # `.moved = True`, driving the gap in Part 2's per-chromosome
        # row), UNLESS one side is trivial (0 pieces), in which case the
        # substantive side always moves (a bare telomere tip has no home
        # to leave).
        bp_u, bp_v = node_u[0], node_v[0]
        M_bps, N_bps = _unit_bps(M), _unit_bps(N)
        is_local = (bp_u in N_bps) or (bp_v in M_bps)
        if not is_local:
            if M.pieces and not N.pieces:
                mover = M
            elif N.pieces and not M.pieces:
                mover = N
            elif M.pieces and N.pieces:
                mover = M if bp_position[bp_u] < bp_position[bp_v] else N
            else:
                mover = None  # both trivial -- nothing to mark
            if mover is not None:
                for p in mover.pieces:
                    p.moved = True

        if M.pieces:
            if M.right_free == node_u:
                pass
            elif M.left_free == node_u:
                M = M.reversed()
            else:
                raise AssertionError(f'node_u {node_u} not a free end of its own molecule (rearr {rearr_num})')
        if N.pieces:
            if N.left_free == node_v:
                pass
            elif N.right_free == node_v:
                N = N.reversed()
            else:
                raise AssertionError(f'node_v {node_v} not a free end of its own molecule (rearr {rearr_num})')

        for end in {M.left_free, M.right_free, N.left_free, N.right_free}:
            free_end_to_molecule.pop(end, None)

        # BUG FIXED 2026-08-13 (caught testing against P05-1657): when M or
        # N is a TRIVIAL (0-piece) telomere-facing tip, its left_free and
        # right_free are BOTH just the tip node itself. Naively using
        # M.left_free / N.right_free as the merged molecule's new free end
        # would re-expose the just-consumed node for future merges (wrong
        # -- it has nothing left to offer) and would also mislabel a
        # resolved-but-chain-external terminus as if it were still open.
        # A consumed trivial tip must become a non-mergeable
        # ('CHAIN_TERMINAL', node) sentinel instead, carrying the node it
        # terminated at (for accurate display, e.g. "...fused to BP19(-),
        # nothing further in this chain's data") -- distinct from a
        # never-touched trivial tip, which keeps its real node as its free
        # end (still open, still eligible for a later merge).
        merged_left = M.left_free if M.pieces else ('CHAIN_TERMINAL', node_u)
        merged_right = N.right_free if N.pieces else ('CHAIN_TERMINAL', node_v)

        merged = Molecule(M.pieces + N.pieces, merged_left, merged_right)
        left_is_terminal = isinstance(merged_left, tuple) and merged_left[0] == 'CHAIN_TERMINAL'
        right_is_terminal = isinstance(merged_right, tuple) and merged_right[0] == 'CHAIN_TERMINAL'
        if not left_is_terminal:
            free_end_to_molecule[merged_left] = merged
        if not right_is_terminal:
            free_end_to_molecule[merged_right] = merged
        # BUG FIXED 2026-08-13 (caught running the full 194-chain cohort,
        # not the P05-1657 test -- P05-1657 never happens to produce this
        # case): when BOTH ends of a merge become CHAIN_TERMINAL (the
        # molecule connects out on both sides with nothing further to
        # attach, e.g. two independent trivial-tip fusions on either end
        # of a small internal piece), NEITHER end gets a live dict entry
        # -- the molecule becomes unreachable and silently vanishes from
        # `open_molecules` at the end, along with every piece it
        # contains. Verified concretely on PR-3042 chain 10 (chr13):
        # pieces (446,409) and (971,633) disappeared entirely, causing a
        # KeyError downstream in build_chromosome_rows. Fixed by keeping
        # an explicit side list for exactly this case.
        if left_is_terminal and right_is_terminal:
            fully_capped_molecules.append(merged)
        rearr_log.append((rearr_num, node_u, node_v, 'MERGED'))

    open_molecules = list(fully_capped_molecules)
    seen_ids = {id(m) for m in fully_capped_molecules}
    for mol in free_end_to_molecule.values():
        if id(mol) not in seen_ids:
            seen_ids.add(id(mol))
            open_molecules.append(mol)

    return {
        'open_molecules': open_molecules,
        'closed_molecules': closed_molecules,
        'rearr_log': rearr_log,
        'n_unresolved_rearr_ends': n_unresolved,
    }


def build_piece_lookup(result):
    """(chrom, low_bp, high_bp) -> Piece, across ALL final molecules
    (open AND closed) -- used by Part 3 to look up a DB-span segment's
    final `.moved`/`.flipped` status directly, without re-deriving the
    chromosome-row layout. Every native reference segment is guaranteed
    to appear exactly once here (the assembly algorithm never destroys a
    piece, only relocates/reorients it -- verified empirically 2026-08-13
    via n_unresolved_rearr_ends == 0 across all 194 DB-bearing chains,
    see DELETION_BRIDGE_ANALYSIS_PLAN.md's Part 3 section)."""
    lookup = {}
    for mol in result['open_molecules']:
        for p in mol.pieces:
            lookup[(p.chrom, p.low_bp, p.high_bp)] = p
    for _, mol in result['closed_molecules']:
        for p in mol.pieces:
            lookup[(p.chrom, p.low_bp, p.high_bp)] = p
    return lookup


def build_chromosome_rows(chain_rows, result):
    """
    Reorganizes the assembled molecules (from assemble_final_structure)
    into ONE ROW PER ORIGINAL CHROMOSOME, matching Part 1's own
    chromosome-track layout and the user's own hand-drawn notation style
    (confirmed 2026-08-13) -- NOT one row per disconnected final
    molecule (an earlier, since-corrected design).

    Walks each chromosome's native breakpoint order. A native segment
    whose piece was never `.moved` becomes an ANCHOR point: the WHOLE
    final molecule it belongs to is rendered there (this naturally pulls
    in any relocated content, e.g. Piece 1 appearing mid-molecule,
    without extra bookkeeping, since the piece/flip computation already
    handles that). A `.moved` segment shows as a GAP at its native slot
    instead -- its real content appears later, wherever its molecule's
    anchor point is (same chromosome, further along the row, OR a
    DIFFERENT chromosome's row if it fused across chromosomes, e.g.
    P05-1657's [3-7] piece reappearing on the chr12 row). Bookend
    telomere-facing tips are treated the same way: genuinely untouched
    ones render as isolated tips; ones consumed by a real fusion
    (identified via the CHAIN_TERMINAL sentinel) become the anchor point
    for whatever molecule attaches there, entered "from the outside in"
    (piece order reversed for display only -- this does NOT alter the
    already-computed, apostrophe-driving `.flipped` value on any piece,
    only which end of the molecule is listed first on this specific row).

    Returns: dict chrom -> ordered list of row items:
      ('gap', {'chrom','low_bp','high_bp'})
      ('molecule', Molecule, is_closed, reversed_for_display)
      ('isolated_tip', node)
    """
    piece_to_mol = {}
    for mol in result['open_molecules']:
        for p in mol.pieces:
            piece_to_mol[(p.chrom, p.low_bp, p.high_bp)] = (mol, False)
    for _, mol in result['closed_molecules']:
        for p in mol.pieces:
            piece_to_mol[(p.chrom, p.low_bp, p.high_bp)] = (mol, True)

    node_to_mol = {}
    chain_terminal_to_mol = {}
    for mol in result['open_molecules']:
        for end in (mol.left_free, mol.right_free):
            if isinstance(end, tuple) and len(end) == 2 and end[0] == 'CHAIN_TERMINAL':
                chain_terminal_to_mol[end[1]] = mol
            else:
                node_to_mol[end] = mol

    rendered = set()  # (chrom, low_bp, high_bp) piece keys already shown

    def _mark_rendered(mol):
        for p in mol.pieces:
            rendered.add((p.chrom, p.low_bp, p.high_bp))

    def _render_tip(node, reversed_if_terminal):
        """Always returns a LIST starting with a ('telomere', node)
        marker -- the real reference stretch leading up to this
        outermost breakpoint from the true chromosome end always exists
        physically regardless of what happened AT the breakpoint itself
        (added 2026-08-13 per user feedback: these telomere-ward
        stretches -- e.g. TELOMERE----BP2(-) -- were previously only a
        text label with no visual bar). If this exact node was also
        consumed by a real fusion, a ('molecule', ...) entry follows it
        (rendered starting from this tip, i.e. entered "from the outside
        in" -- reversed_if_terminal orients a CHAIN_TERMINAL-anchored
        molecule accordingly); if genuinely untouched, the telomere
        marker is the whole entry."""
        if node in node_to_mol:
            mol = node_to_mol[node]
            if not mol.pieces:
                return [('telomere', node, False)]
            _mark_rendered(mol)
            return [('telomere', node, True), ('molecule', mol, False, False)]
        if node in chain_terminal_to_mol:
            mol = chain_terminal_to_mol[node]
            _mark_rendered(mol)
            return [('telomere', node, True), ('molecule', mol, False, reversed_if_terminal)]
        return [('telomere', node, False)]  # defensive fallback; should not normally trigger

    rows = {}
    for chrom, grp in chain_rows.groupby('chromosome'):
        chrom = int(chrom)
        bps_sorted = [int(bp) for bp in
                      grp[['Breakpoint number', 'position']]
                      .drop_duplicates(subset=['Breakpoint number'])
                      .sort_values(['position', 'Breakpoint number'])['Breakpoint number']]
        row = list(_render_tip((bps_sorted[0], '-'), reversed_if_terminal=True))
        for i in range(len(bps_sorted) - 1):
            low_bp, high_bp = bps_sorted[i], bps_sorted[i + 1]
            key = (chrom, low_bp, high_bp)
            if key in rendered:
                continue
            mol, is_closed = piece_to_mol[key]
            piece_obj = next(p for p in mol.pieces if (p.chrom, p.low_bp, p.high_bp) == key)
            if piece_obj.moved:
                row.append(('gap', {'chrom': chrom, 'low_bp': low_bp, 'high_bp': high_bp}))
            else:
                _mark_rendered(mol)
                row.append(('molecule', mol, is_closed, False))
        row.extend(_render_tip((bps_sorted[-1], '+'), reversed_if_terminal=False))
        rows[chrom] = row
    return rows


def _end_label(end, is_left):
    """Render one end of a molecule: a real, still-open (bp,strand) node
    shows as 'TELOMERE' (matches the project-wide convention that 'no
    further recorded breakpoint in this chain' = telomere sentinel,
    applied here identically); a ('CHAIN_TERMINAL', node) sentinel shows
    the specific node it fused to before running out of chain data, e.g.
    a cross-chromosome jump to a breakpoint with nothing further known
    beyond it (matches the user's own P05-1657 notation, e.g.
    'TELOMERE----19[-] 7[-]----3[+]' — BP19 is shown explicitly, telomere
    only prefixes it)."""
    if isinstance(end, tuple) and len(end) == 2 and end[0] == 'CHAIN_TERMINAL':
        bp, strand = end[1]
        node_txt = f'{bp}[{strand}]'
        return f'TELOMERE----{node_txt}' if is_left else f'{node_txt}----TELOMERE'
    return 'TELOMERE'


def molecule_to_text(mol):
    """Human-readable left-to-right rendering, e.g.
    'TELOMERE----6[+]----18[-]  18[+]----15[-]  3[-]----2[+]  14[-]----15[+]----TELOMERE'
    matching the user's own by-hand notation style (two-space gap = real
    fusion junction, ---- = intact reference segment)."""
    if not mol.pieces:
        # A lone, never-consumed telomere-facing tip -- a single node, not
        # a segment (nothing on its "other side" since it's zero-length).
        bp, strand = mol.left_free
        return f'TELOMERE----{bp}[{strand}]  (isolated -- no rearrangement or reference partner in this chain)'
    parts = [p.current_span_text() for p in mol.pieces]
    left = _end_label(mol.left_free, is_left=True)    # 'TELOMERE' or 'TELOMERE----19[-]'
    right = _end_label(mol.right_free, is_left=False)  # 'TELOMERE' or '19[-]----TELOMERE'
    return left + '  ' + '  '.join(parts) + '  ' + right

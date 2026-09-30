"""Independent replay of a cascade certificate (certificate.json).

usage: ~/.venvs/atlas/bin/python cascade_check.py <certificate.json>
           [--farsidediagram PATH] [--knot-table CSV] [--link-table CSV]

Written separately from the C++ cascade (cascade/*.cpp), in the spirit of
frontier.py --check: nothing here calls the cascade's composition, identity
or bookkeeping code. What it does use:
  - farsidediagram --gauss (the search pipeline's own far-side drawer,
    validated separately: knotbuilder/README.md) to re-read each witness:
    from its pair signature, or -- an in-process hop's witness -- from its
    triangles in the row's thickening (--faces), which farsidediagram
    rebuilds with the search's own embedding checks, after the rebuilt
    thickening's digest has been matched to the one the certificate records;
  - Regina (Link.fromData / fromPD) and SnapPy (isometries) as libraries.

For every record of the proof it re-derives the claim:
  leaf            the unknot really is one; a literature value is the table's
                  (and not the target's own); a direct witness bounds the row
                  alone with that grouping and genus
  witness-*       the witness's surface components on both ends (from
                  farsidediagram), the row -> node component map (its own
                  diagram-isomorphism search), the far side's split pieces
                  (its own union-find on crossings) and each piece's identity
                  with its node (a diagram isomorphism under the certificate's
                  component map, after its own simplification and removal of
                  nugatory crossings, or a SnapPy isometry carrying meridians to
                  meridians with one sign realising that map); then the
                  glued surface's partition and genus by an Euler-
                  characteristic count
  split-*         pieces combine in disjoint balls / restrict, with the maps
                  checked against the far side's pieces

Prints one line per record and a verdict. Exit 0 iff every record checks.
"""
import argparse, csv, hashlib, json, os, subprocess, sys
from collections import defaultdict
from itertools import permutations

csv.field_size_limit(10**9)
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_FSD = os.path.normpath(os.path.join(
    HERE, '../../../../build/utils/surfer/farsidediagram'))
ATLAS = os.environ.get('CASCADE_ATLAS', '/home/john/Projects/triangles/cobordism-atlas')

import regina  # noqa: E402

# ---------------------------------------------------------------- diagrams


class Gauss:
    """Signed Gauss data: signs[k] of crossing k; comps[c] = [+-(k+1), ...]
    (+ over, - under), in each component's order."""

    def __init__(self, signs, comps):
        self.signs = list(signs)
        self.comps = [list(c) for c in comps]

    @staticmethod
    def of_link(link):
        signs = [link.crossing(k).sign() for k in range(link.size())]
        comps = []
        for i in range(link.countComponents()):
            s = link.component(i)
            word = []
            if s:
                start = s
                while True:
                    k = s.crossing().index() + 1
                    word.append(k if s.strand() == 1 else -k)
                    s = s.next()
                    if s == start:
                        break
            comps.append(word)
        return Gauss(signs, comps)

    def link(self):
        return regina.Link.fromData(self.signs, self.comps) if self.comps else regina.Link()

    def lk(self):
        m = len(self.comps)
        owner = {}
        for c, w in enumerate(self.comps):
            for x in w:
                owner[(abs(x) - 1, x > 0)] = c
        out = [[0] * m for _ in range(m)]
        for k, s in enumerate(self.signs):
            a, b = owner[(k, True)], owner[(k, False)]
            if a != b:
                out[a][b] += s
                out[b][a] += s
        return [[v // 2 for v in row] for row in out]

    def mirror(self):
        return Gauss([-s for s in self.signs], [[-x for x in w] for w in self.comps])

    def reverse_all(self):
        return Gauss(self.signs, [list(reversed(w)) for w in self.comps])

    def turned_over(self):
        """The whole diagram turned over (a half turn about an axis in the
        plane): every crossing swaps over and under, every sign stays."""
        return Gauss(self.signs, [[-x for x in w] for w in self.comps])


def nugatory_side(g, c):
    """For a self-crossing c: the crossings cut off with the loop its
    component runs between its two visits to c, if deleting c disconnects
    that loop (and whatever crosses it) from the component's other loop;
    else None. By union-find over the projection graph with c deleted."""
    where = [(ci, p) for ci, w in enumerate(g.comps) for p, x in enumerate(w) if abs(x) - 1 == c]
    if len(where) != 2 or where[0][0] != where[1][0]:
        return None
    w = g.comps[where[0][0]]
    i, j = where[0][1], where[1][1]
    parent = list(range(len(g.signs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for word in g.comps:
        for p in range(len(word)):
            a, b = abs(word[p]) - 1, abs(word[(p + 1) % len(word)]) - 1
            if c not in (a, b):
                parent[find(a)] = find(b)
    inner = {abs(x) - 1 for x in w[i + 1:j]}
    outer = {abs(x) - 1 for x in w[j + 1:] + w[:i]}
    roots = {find(k) for k in inner}
    if roots & {find(k) for k in outer}:
        return None
    return {k for k in range(len(g.signs)) if k != c and find(k) in roots}


def lift_split(g):
    """g with every component that passes over (or under) at every crossing
    it meets lifted off, leaving it crossingless: having no self-crossing, it
    is an unknot above (below) the rest, so lifting it is an isotopy that
    splits it off. Repeated until none is left."""
    while True:
        lift = next((c for c, w in enumerate(g.comps)
                     if w and (all(x > 0 for x in w) or all(x < 0 for x in w))), None)
        if lift is None:
            return g
        gone = {abs(x) - 1 for x in g.comps[lift]}
        keep = [k for k in range(len(g.signs)) if k not in gone]
        new = {k: i + 1 for i, k in enumerate(keep)}
        g = Gauss([g.signs[k] for k in keep],
                  [[(1 if x > 0 else -1) * new[abs(x) - 1] for x in w if abs(x) - 1 not in gone]
                   for w in g.comps])


def remove_nugatory(g):
    """g with every nugatory crossing removed, each by turning over the side
    it cuts off (a half turn about an axis in the plane through it, a rigid
    motion): the crossing goes, that side's crossings swap over and under,
    and every sign stays. An isotopy that keeps component order and
    orientation."""
    while True:
        c = next((k for k in range(len(g.signs)) if nugatory_side(g, k) is not None), None)
        if c is None:
            return g
        side = nugatory_side(g, c)
        comps = []
        for w in g.comps:
            nw = []
            for x in w:
                k = abs(x) - 1
                if k == c:
                    continue
                label = (k if k < c else k - 1) + 1
                nw.append(label if (x > 0) != (k in side) else -label)
            comps.append(nw)
        g = Gauss([s for k, s in enumerate(g.signs) if k != c], comps)


def pd_to_link(pd_text):
    nums, cur = [], ''
    for ch in pd_text:
        if ch.isdigit():
            cur += ch
        elif cur:
            nums.append(int(cur)); cur = ''
    if cur:
        nums.append(int(cur))
    shift = 1 if min(nums) == 0 else 0
    xs = [[nums[i] + shift, nums[i + 1] + shift, nums[i + 2] + shift, nums[i + 3] + shift]
          for i in range(0, len(nums), 4)]
    return regina.Link.fromPD(xs)


def iso_with_map(a, b, comp_map):
    """Is there a diagram isomorphism a -> b (orientation kept, no mirror)
    sending component c of a to comp_map[c] of b? Exhaustive over each
    component's starting point."""
    if len(a.comps) != len(b.comps) or len(a.signs) != len(b.signs):
        return False
    cross = {}

    def rec(c):
        if c == len(a.comps):
            return True
        wa, wb = a.comps[c], b.comps[comp_map[c]]
        if len(wa) != len(wb):
            return False
        if not wa:
            return rec(c + 1)
        for rot in range(len(wa)):
            added, ok = [], True
            for i, xa in enumerate(wa):
                xb = wb[(i + rot) % len(wb)]
                if (xa > 0) != (xb > 0) or a.signs[abs(xa) - 1] != b.signs[abs(xb) - 1]:
                    ok = False; break
                ka, kb = abs(xa) - 1, abs(xb) - 1
                if ka in cross:
                    if cross[ka] != kb:
                        ok = False; break
                else:
                    if kb in cross.values():
                        ok = False; break
                    cross[ka] = kb; added.append(ka)
            if ok and rec(c + 1):
                return True
            for ka in added:
                del cross[ka]
        return False

    return rec(0)


def iso_maps(a, b):
    """Every component map of a diagram isomorphism a -> b (orientation kept,
    no mirror): one per component permutation some isomorphism realises."""
    m = len(a.comps)
    if m > 7:
        raise RuntimeError('too many components for an exhaustive map search')
    for perm in permutations(range(m)):
        if iso_with_map(a, b, list(perm)):
            yield list(perm)


def find_iso(a, b, allow_mirror, allow_reverse):
    """Some (map, mirrored, reversed) with iso_with_map, trying every map."""
    m = len(a.comps)
    if m > 7:
        raise RuntimeError('too many components for an exhaustive map search')
    for mir in ([False, True] if allow_mirror else [False]):
        for rev in ([False, True] if allow_reverse else [False]):
            t = a.mirror() if mir else a
            t = t.reverse_all() if rev else t
            for perm in permutations(range(m)):
                if iso_with_map(t, b, list(perm)):
                    return list(perm), mir, rev
    return None


def split_pieces(g):
    """Components joined by crossings, as lists of component indices."""
    m = len(g.comps)
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    owner = defaultdict(list)
    for c, w in enumerate(g.comps):
        for x in w:
            owner[abs(x) - 1].append(c)
    for k, cs in owner.items():
        for c in cs[1:]:
            parent[find(c)] = find(cs[0])
    groups = defaultdict(list)
    for c in range(m):
        groups[find(c)].append(c)
    return [sorted(v) for v in groups.values()]


def sub_gauss(g, comps):
    """The sub-diagram of the given components (a split piece), crossings
    renumbered; component i of the result is comps[i]."""
    keep = sorted({abs(x) - 1 for c in comps for x in g.comps[c]})
    new = {k: i for i, k in enumerate(keep)}
    return Gauss([g.signs[k] for k in keep],
                 [[(1 if x > 0 else -1) * (new[abs(x) - 1] + 1) for x in g.comps[c]] for c in comps])


def is_unknot_piece(g):
    if len(g.comps) != 1:
        return False
    if not g.comps[0]:
        return True
    l = g.link()
    l.simplify()
    if l.size() == 0:
        return True
    try:
        return l.isTrivial()
    except Exception:
        return False


# The atlas's own Python identification code (a separate implementation of
# what the C++ does, validated against it): tools/row_certificates.py and
# tools/farside/diagram/fsid.py. Imported lazily; nothing here calls their
# main().
_atlas_mods = {}


def atlas_module(name):
    if name not in _atlas_mods:
        sys.path.insert(0, os.path.join(ATLAS, 'tools'))
        sys.path.insert(0, os.path.join(ATLAS, 'tools', 'farside', 'diagram'))
        _atlas_mods[name] = __import__(name)
    return _atlas_mods[name]


def snappy_link(g):
    """SnapPy's Link for a Gauss diagram, components in the same order, each
    oriented as in the diagram (via Regina's PD, labels along components)."""
    import snappy
    l = g.link()
    return snappy.Link([list(x) for x in l.pdData()])


def same_link_by_isometry(piece, node, comp_map, mirrored, reversed_):
    """An isometry of the complements carrying meridians to meridians with
    ONE sign on every component (row_certificates.meridian_signs, the
    atlas's check), realising comp_map (piece component i -> node component
    comp_map[i]). Returns (ok, why)."""
    import snappy
    rc = atlas_module('row_certificates')
    try:
        # Cusps in COMPONENT order, meridians the diagram's own (Regina's
        # SnapPea construction, as exactnaming uses it): the cusp permutation
        # is then comparable with comp_map. SnapPy's Link may reorder
        # components, so it is not used here.
        mp = snappy.Manifold(regina.SnapPeaTriangulation(piece.link()).snapPea())
        mn = snappy.Manifold(regina.SnapPeaTriangulation(node.link()).snapPea())
    except Exception as e:
        return False, f'no complement: {e}'
    for attempt in range(4):
        try:
            isos = mp.is_isometric_to(mn, return_isometries=True)
            break
        except Exception:
            mp.randomize()
    else:
        return False, 'SnapPy could not decide (not hyperbolic?)'
    for iso in isos:
        if list(iso.cusp_images()) != list(comp_map):
            continue
        signs = rc.meridian_signs(iso)
        if signs and len(set(signs)) == 1:
            return True, 'isometry'
    return False, f'no meridian-preserving isometry realising the map ({len(isos)} isometries)'


def prove_table_identity(g, name, tables):
    """Re-prove that diagram g IS the table entry `name` (as an oriented link
    up to mirror and global reversal), with the atlas's code. Knots: fsid's
    exterior isometry against OUR table PD (Gordon-Luecke). Links: an
    isometry carrying meridians with one sign onto the entry's own table PD
    (row_certificates.certify_hyperbolic); or, when every orientation variant
    of the base has the same literature value, fsid's base identification."""
    if len(g.comps) == 1:
        fsid = atlas_module('fsid')
        got, proof, _ = fsid.identify_knot(snappy_link(g))
        base = name.lstrip('m')
        return (got is not None and got.lstrip('m') == base), f'{proof}: {got}'
    import snappy
    rc = atlas_module('row_certificates')
    pd = tables['pd'].get(name)
    if pd is None:
        return False, f'{name} has no table PD'
    ref = snappy.Link(rc.pd_list(pd)).exterior()
    ours = snappy_link(g).exterior()
    if rc.hyperbolic(ours):
        cert = rc.certify_hyperbolic(ours, ref)
        if cert and cert.get('sign') is not None:
            return True, 'isometry with uniform meridian signs'
        return False, f'no uniform-sign meridian isometry: {cert}'
    base = name.split('{')[0]
    variants = {v for n, v in tables['g4'].items() if n.split('{')[0] == base}
    if len(variants) == 1:
        fsid = atlas_module('fsid')
        got, proof, _ = fsid.identify_link(snappy_link(g))
        return got == base, f'{proof}: {got} (every variant of {base} has g4 {variants})'
    return False, 'non-hyperbolic link whose variants differ: orientation not provable here'


# ---------------------------------------------------------------- gluing


def glue(shape, glued_side, glued_labels, glued_genus):
    """Independent Euler-characteristic model of capping one side of a
    cobordism: returns (labels of the free side's curves, genus upper bound).
    shape: {'components', 'genus', 'inComponent', 'outComponent'}."""
    ins, outs = shape['inComponent'], shape['outComponent']
    glued_c, free_c = (outs, ins) if glued_side == 'out' else (ins, outs)
    nc = shape['components']
    blocks = sorted(set(glued_labels))
    bidx = {b: i for i, b in enumerate(blocks)}
    nv = nc + len(blocks)
    adj = defaultdict(list)
    for j, c in enumerate(glued_c):
        u, v = c, nc + bidx[glued_labels[j]]
        adj[u].append(v); adj[v].append(u)
    comp = [-1] * nv
    ncomp = 0
    for s in range(nv):
        if comp[s] >= 0:
            continue
        stack = [s]; comp[s] = ncomp
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if comp[y] < 0:
                    comp[y] = ncomp; stack.append(y)
        ncomp += 1
    # circles per vertex, and chi per graph component, with genus spread as
    # all of C's on vertex 0 and all of G's on the first block (the total is
    # what matters: chi is additive)
    circles = [0] * nv
    for c in ins + outs:
        circles[c] += 1
    for lab in glued_labels:
        circles[nc + bidx[lab]] += 1
    genus_v = [0] * nv
    genus_v[0] += shape['genus']
    if blocks:
        genus_v[nc] += glued_genus
    chi = [0] * ncomp
    for v in range(nv):
        chi[comp[v]] += 2 - 2 * genus_v[v] - circles[v]
    free_count = [0] * ncomp
    for c in free_c:
        free_count[comp[c]] += 1
    total = 0
    for k in range(ncomp):
        twice = 2 - free_count[k] - chi[k]
        assert twice % 2 == 0 and twice >= 0, 'non-integral genus'
        total += twice // 2
    return [comp[c] for c in free_c], total


def normalize(labels):
    seen = {}
    return [seen.setdefault(x, len(seen)) for x in labels]


def refines(fine, coarse):
    """Whether every block of `fine` lies inside a block of `coarse`."""
    block_of = {}
    for i, b in enumerate(fine):
        if block_of.setdefault(b, coarse[i]) != coarse[i]:
            return False
    return True


def partition_labels(text):
    """'{0,2}{1}' -> labels per element."""
    blocks = [b for b in text.strip('{}').split('}{')] if text != '{}' else []
    lab = {}
    for i, b in enumerate(blocks):
        for x in (b.split(',') if b else []):
            lab[int(x)] = i
    return [lab[i] for i in range(len(lab))]


# ---------------------------------------------------------------- replay

class Checker:
    def __init__(self, cert, fsd, tables):
        self.cert = cert
        self.fsd = fsd
        self.tables = tables
        self.records = {r['id']: r for r in cert['records']}
        self.nodes = {n['id']: n for n in cert['nodes']}
        self.fsd_cache = {}
        self.split_claims = {}  # split whole node -> the pieces a witness drew
        self.problems = []
        self.notes = []

    def is_target_link(self, name):
        """Whether the table entry `name` is the target's own link (one
        oriented link up to mirror and global reversal) under another name:
        by the class table when there is one, else by the meridian-carrying
        isometry that proves every identity here (one-sided: an identity it
        cannot prove is not refused). A literature leaf on such an entry would
        prove the target by its own literature value."""
        target = self.cert['target']
        classes = self.tables.get('classes', {})
        if classes.get(name, name) == classes.get(target, target):
            return True, 'data/table_link_classes.csv'
        if name not in self.tables['pd'] or target not in self.tables['pd']:
            return False, ''
        tnode = [n['id'] for n in self.cert['nodes'] if n['label'].startswith('target')][0]
        try:
            return prove_table_identity(self.node_gauss(tnode), name, self.tables)
        except Exception as e:  # noqa: BLE001 -- unproved is not refused
            return False, f'{type(e).__name__}: {e}'

    def node_gauss(self, n):
        node = self.nodes[n]
        if 'gauss' in node:
            # The node's own diagram: component maps refer to its order.
            return Gauss(node['signs'], node['gauss'])
        if 'pd' not in node:
            return Gauss([], [[]] * node['components'])
        self.notes.append(f'node {n}: no Gauss data; component order taken from its PD')
        return Gauss.of_link(pd_to_link(node['pd']))

    def pairsig(self, hop_dir, key):
        with open(os.path.join(hop_dir, 'cob.csv')) as fh:
            for r in csv.DictReader(fh):
                if hashlib.sha1(r['pairsig'].encode()).hexdigest()[:12] == key:
                    return r
        raise RuntimeError(f'witness {key} not in {hop_dir}/cob.csv')

    def redraw(self, hop_dir, row_pd, key, record=None):
        ck = (hop_dir, key)
        if ck in self.fsd_cache:
            return self.fsd_cache[ck]
        faces = record is not None and 'faces' in record
        if faces:
            # An in-process witness: its triangles in the row's thickening,
            # rebuilt by farsidediagram --faces with the search's own checks.
            # The digest below proves the rebuilt thickening is the one they
            # index.
            w = {'data': ','.join(str(t) for t in record['faces']),
                 'thicken_layers': str(record.get('layers', 2)), 'genus': None}
        elif record is not None and 'pairsig' in record:
            # A master witness: its pair signature travels in the certificate,
            # under a provenance-prefixed key ("master:<sha1[:12]>").
            ps = record['pairsig']
            if hashlib.sha1(ps.encode()).hexdigest()[:12] != key.split(':')[-1]:
                raise RuntimeError(f'inline pair signature does not hash to {key}')
            w = {'data': ps, 'thicken_layers': str(record.get('layers', 2)),
                 'genus': None}
        else:
            w = self.pairsig(hop_dir, key)
            w['data'] = w['pairsig']
        out = subprocess.run([self.fsd, '--layers', w['thicken_layers'] or '2', '--gauss']
                             + (['--faces'] if faces else []) + [row_pd],
                             input=f'0 {w["data"]}\n', capture_output=True,
                             text=True).stdout.splitlines()
        row = dict(kv.split('=', 1) for kv in out[0].split(' ')[1:])
        if faces and row.get('build') != record.get('build'):
            raise RuntimeError(f"the rebuilt thickening's digest {row.get('build')} is not "
                               f"the one its faces were recorded in ({record.get('build')})")
        wl = [l for l in out if l.startswith('W ')][0]
        if ' ok ' not in wl:
            raise RuntimeError(f'farsidediagram: {wl}')
        f = dict(kv.split('=', 1) for kv in wl.split(' ', 3)[3].split(' '))
        if 'genus' not in f:
            raise RuntimeError('farsidediagram gave no genus (rebuild it)')
        surface_genus = int(f['genus'])  # from the pair signature's own surface
        if w['genus'] is not None and int(w['genus']) != surface_genus:
            raise RuntimeError(f'witness file genus {w["genus"]} != surface genus {surface_genus}')
        res = {
            'genus': surface_genus,
            'row': Gauss(json.loads(row['signs']), json.loads(row['gauss'])),
            'far': Gauss(json.loads(f['signs']), json.loads(f['gauss'])),
            'surface': json.loads(f['surface']),
            'incoming': {int(k): v for k, v in json.loads(f['incoming']).items()},
            'edges': json.loads(f['edges']) if 'edges' in f else None,
        }
        self.fsd_cache[ck] = res
        return res

    def same_as_node(self, sub, node, comp_map, mirrored, reversed_):
        """Is diagram `sub` the node's link under comp_map (and flags)?
        A diagram match of `sub` or one of our own simplifications of it,
        with and without its nugatory crossings removed, either way up
        (each an isotopy, keeping component order); else an isometry
        carrying meridians with one sign realising comp_map; for a knot,
        any exterior isometry (Gordon-Luecke) or fsid proving both the same
        table knot."""
        for attempt in range(41):
            cand = sub
            if attempt:
                l = sub.link()
                l.simplify()
                cand = Gauss.of_link(l)
            reduced = remove_nugatory(cand)
            for c in (cand, reduced, reduced.turned_over()):
                t = c.mirror() if mirrored else c
                t = t.reverse_all() if reversed_ else t
                if iso_with_map(t, node, comp_map):
                    return True, 'diagram'
        ok, why = same_link_by_isometry(sub, node, comp_map, mirrored, reversed_)
        if ok or len(sub.comps) != 1:
            return ok, why
        import snappy
        try:
            a, b = snappy_link(sub).exterior(), snappy_link(node).exterior()
            if a.is_isometric_to(b):
                return True, 'knot exterior isometry'
        except Exception:
            pass
        fsid = atlas_module('fsid')
        na, pa, _ = fsid.identify_knot(snappy_link(sub))
        nb, pb, _ = fsid.identify_knot(snappy_link(node))
        if na is not None and na.lstrip('m') == (nb or '').lstrip('m'):
            return True, f'both proved {na} ({pa}/{pb})'
        return False, f'{why}; knot not identified ({na}, {nb})'

    def fail(self, rid, why):
        self.problems.append(f'record {rid}: {why}')
        return False

    # -- leaves
    def check_leaf(self, r):
        src = r['source']
        node = self.nodes[r['node']]
        if src.startswith('unknot'):
            g = self.node_gauss(r['node'])
            if node['components'] != 1 or not is_unknot_piece(g):
                return self.fail(r['id'], 'unknot leaf on a node that is not the unknot')
            return r['genus'] >= 0 and r['partition'] == '{0}'
        if src.startswith('literature '):
            _, name, g4 = src.split(' ', 2)
            t = self.tables['g4'].get(name)
            if t is None:
                return self.fail(r['id'], f'{name} not in the tables')
            hi = int(t.strip('[]').split(';')[-1])
            if hi != r['genus'] or t != g4:
                return self.fail(r['id'], f'table says {t}, record says {g4}/{r["genus"]}')
            if name == self.cert['target']:
                return self.fail(r['id'], "the target's own literature value")
            own, why_own = self.is_target_link(name)
            if own:
                return self.fail(r['id'], f"{name} is the target's own link ({why_own}), "
                                          "so its literature value is the target's own")
            if len(set(partition_labels(r['partition']))) != 1:
                return self.fail(r['id'], 'literature bounds only the connected genus')
            ok, why = prove_table_identity(self.node_gauss(r['node']), name, self.tables)
            if not ok:
                return self.fail(r['id'], f'node {r["node"]} not proved to be {name}: {why}')
            self.notes.append(f'record {r["id"]}: literature {name} {t}; identity: {why}')
            return True
        if src.startswith('direct witness '):
            d = self.redraw(r['hop_dir'], r['row_pd'], r['witness'], r)
            if d['far'].comps:
                return self.fail(r['id'], 'a direct witness with a far side')
            return self.fail(r['id'], 'direct witnesses: replay not implemented')
        return self.fail(r['id'], f'unknown leaf source {src!r}')

    # -- witnesses
    def check_witness(self, r):
        return self.check_witness_edge(r) and self.check_glue(r)

    def check_witness_edge(self, r):
        """The edge alone: the surface redrawn, the row -> node map, the
        shape, the far side's pieces and their identities. Shared by an
        upper record (which then glues) and a lower fact (which then
        transports)."""
        d = dict(self.redraw(r['hop_dir'], r['row_pd'], r['witness'], r))
        rid = r['id']
        # Put our redraw's far-side curves in the CERTIFICATE's order: two
        # reads may list them differently; a curve is its set of edges.
        if 'farCurveEdges' in r:
            if d['edges'] is None:
                return self.fail(rid, 'farsidediagram gave no edges (rebuild it)')
            where = {tuple(e): j for j, e in enumerate(d['edges'])}
            perm = [where.get(tuple(e)) for e in r['farCurveEdges']]
            if None in perm or sorted(perm) != list(range(len(d['edges']))):
                return self.fail(rid, 'far-side curves do not match by edge set')
            d['surface'] = [d['surface'][perm[j]] for j in range(len(perm))]
            d['far'] = Gauss(d['far'].signs, [d['far'].comps[perm[j]] for j in range(len(perm))])
        # the row -> node map. Step 1: the redrawn row IS the row's own
        # diagram (its PD, read by Regina), by our own isomorphism search.
        # Step 2: that diagram IS the node, under the certificate's
        # row_node_map (identity when the row is the node's own diagram),
        # proved like a far-side piece.
        in_node = r['in']
        node_g = self.node_gauss(in_node)
        # When the row's diagram has automorphisms there are several such
        # maps. Each comes from a symmetry h of (S^3, L); h x id carries the
        # witness to one whose incoming components are permuted by it, and
        # since h is isotopic to the identity in S^3 the far-side curves go
        # back to themselves: the same surface witnesses the shape under
        # every one of them. So the certificate's shape stands if it is the
        # redraw's under SOME map (the far side is compared exactly).
        if 'row_node_map' in r:
            row_g = Gauss.of_link(pd_to_link(r['row_pd']))
            rmap = r['row_node_map']
            ok, why = self.same_as_node(row_g, node_g, rmap, False, False)
            if not ok:
                return self.fail(rid, f'the row diagram is not node {in_node} under its map: {why}')
            maps = [[rmap[m[i]] for i in range(len(m))] for m in iso_maps(d['row'], row_g)]
            if not maps:
                return self.fail(rid, 'the row does not redraw as its own PD (no isomorphism)')
        else:
            maps = list(iso_maps(d['row'], node_g))
            if not maps:
                return self.fail(rid, 'the row does not redraw as its node (no isomorphism)')
        n_in = len(d['row'].comps)
        # the shape, from the redraw, under each map
        sc_of_row = {}
        for sc, comps in d['incoming'].items():
            for rc in comps:
                sc_of_row[rc] = sc
        cs = r['shape']
        theirs = normalize(cs['inComponent'] + cs['outComponent'])
        mine = None
        for row_to_node in maps:
            in_comp = [None] * n_in
            for rc in range(n_in):
                in_comp[row_to_node[rc]] = sc_of_row[rc]
            mine = normalize(in_comp + d['surface'])
            if mine == theirs:
                break
        # our in/out orders: in = node order; out = drawn curve order (outMap
        # composes them with the far node, checked below)
        if mine != theirs:
            return self.fail(rid, f'shape differs: redrawn {mine}, certificate {theirs}')
        if cs['genus'] != d['genus']:
            return self.fail(rid, f'genus differs: {d["genus"]} vs {cs["genus"]}')
        # pieces and identities
        far = d['far']
        cert_pieces = r.get('pieces', [])
        claimed = sorted(sorted(p['origins']) for p in cert_pieces)
        # A drawn piece may come apart only after Reidemeister moves, or once
        # a component lying above (below) everything is lifted off: find, for
        # every drawn piece, a simplification (our own Regina calls and
        # lifting; each an isotopy) that splits it exactly as claimed.
        diagram_of = {}  # tuple(origins) -> Gauss of that claimed piece
        for drawn in split_pieces(far):
            want = [c for c in claimed if set(c) <= set(drawn)]
            if sorted(x for c in want for x in c) != sorted(drawn):
                return self.fail(rid, f'pieces differ: drawn {drawn}, claimed {claimed}')
            sub = sub_gauss(far, drawn)
            if len(want) == 1:
                diagram_of[tuple(want[0])] = sub
                continue
            found = None
            for attempt in range(41):
                if attempt:
                    l = sub.link()
                    l.simplify()
                    g = Gauss.of_link(l)  # components keep their order (Link.simplify)
                else:
                    g = sub
                g = lift_split(g)
                groups = [sorted(drawn[i] for i in grp) for grp in split_pieces(g)]
                if sorted(groups) == sorted(want):
                    found = g
                    break
            if found is None:
                return self.fail(rid, f'drawn piece {drawn} did not split as {want} '
                                      f'under 40 simplifications')
            for grp in split_pieces(found):
                diagram_of[tuple(sorted(drawn[i] for i in grp))] = sub_gauss(found, grp)
        for p in cert_pieces:
            sub = diagram_of[tuple(sorted(p['origins']))]
            # sub's components are in increasing curve order; the certificate
            # lists origins in its own order: reorder to match.
            order = sorted(p['origins'])
            sub = Gauss(sub.signs, [sub.comps[order.index(o)] for o in p['origins']])
            if p['method'] == 'unknot':
                if not is_unknot_piece(sub):
                    return self.fail(rid, f'piece {p["origins"]} is not an unknot')
                continue
            node = self.node_gauss(p['node'])
            ok, why = self.same_as_node(sub, node, p['componentMap'], p['mirrored'],
                                        p['reversed'])
            if not ok:
                return self.fail(rid, f'piece {p["origins"]} -> node {p["node"]}: {why}')
        # The maps must agree with the pieces: one piece -> outMap is its
        # component map; several -> the far node is a split whole whose split
        # record's pieceMap is the pieces' maps (checked in check_split via
        # self.split_claims).
        if len(cert_pieces) == 1:
            p = cert_pieces[0]
            if r['out' if r['kind'] == 'witness-forward' else 'out'] != p['node']:
                return self.fail(rid, 'the far node is not its one piece')
            for i, o in enumerate(p['origins']):
                if r['outMap'][o] != p['componentMap'][i]:
                    return self.fail(rid, 'outMap disagrees with the piece map')
        else:
            if r['outMap'] != list(range(len(r['outMap']))):
                return self.fail(rid, 'a split far side must map curves to itself')
            self.split_claims[r['out']] = cert_pieces
        return True

    def check_glue(self, r):
        child = self.records[r['children'][0]]
        cs = r['shape']
        fwd = r['kind'] == 'witness-forward'
        glued_map, free_map = (r['outMap'], r['inMap']) if fwd else (r['inMap'], r['outMap'])
        child_labels = partition_labels(child['partition'])
        glued_labels = [child_labels[glued_map[j]] for j in range(len(glued_map))]
        free_labels, genus = glue(cs, 'out' if fwd else 'in', glued_labels, child['genus'])
        node_labels = [None] * len(free_map)
        for i, v in enumerate(free_map):
            node_labels[v] = free_labels[i]
        if normalize(node_labels) != normalize(partition_labels(r['partition'])) or \
                genus != r['genus']:
            return self.fail(r['id'], f'gluing gives {normalize(node_labels)} g{genus}, '
                                      f'record {r["partition"]} g{r["genus"]}')
        return True

    def check_split(self, r):
        # The split edge must be the decomposition some checked witness drew.
        claim = self.split_claims.get(r['whole'])
        if claim is None:
            # the witness creating this whole comes later in id order only if
            # the proof is malformed; look it up directly
            for w in self.records.values():
                if w['kind'].startswith('witness') and w.get('out') == r['whole']:
                    try:
                        self.check_witness(w)
                    except Exception:
                        pass
            claim = self.split_claims.get(r['whole'])
        if claim is None:
            return self.fail(r['id'], 'no witness drew this split far side')
        if [p['node'] for p in claim] != r['pieces']:
            return self.fail(r['id'], 'split pieces differ from the witness pieces')
        for k, p in enumerate(claim):
            for i, o in enumerate(p['origins']):
                if r['pieceMap'][k][p['componentMap'][i]] != o:
                    return self.fail(r['id'], 'pieceMap disagrees with the piece origins')
        if r['kind'] == 'split-combine':
            labels = {}
            off, genus = 0, 0
            for k, cid in enumerate(r['children']):
                c = self.records[cid]
                if c['node'] != r['pieces'][k]:
                    return self.fail(r['id'], 'child on the wrong piece')
                cl = partition_labels(c['partition'])
                for comp, whole in enumerate(r['pieceMap'][k]):
                    labels[whole] = off + cl[comp]
                off += max(cl) + 1 if cl else 0
                genus += c['genus']
            lab = [labels[i] for i in range(len(labels))]
            if normalize(lab) != normalize(partition_labels(r['partition'])) or genus != r['genus']:
                return self.fail(r['id'], 'split combine does not reproduce')
            return True
        return self.fail(r['id'], 'split-restrict: replay not implemented')

    # -- lower bounds (lower_certificate.json: facts, children before parents)
    #
    # A fact says lower(node, partition) >= value: every surface for the
    # node whose partition refines it has at least that genus. Its reason is
    # replayed here with the same tools as an upper record, plus the
    # transport: capping a surface of the fact's partition (genus 0, since
    # only the addition matters) onto the witness must give the other end's
    # recorded partition -- or one refining what its fact is stored for --
    # and value <= from - addition. Bounds only ever rise, so `from` may be
    # below the fact read today; both are checked. Values may be "inf"
    # (no surface has that partition).
    def value_of(self, v):
        return float('inf') if v == 'inf' else int(v)

    def check_lower_fact(self, f):
        fid = f['id']
        kind = f['kind']
        value = self.value_of(f['value'])
        node = self.nodes[f['node']]
        labels = partition_labels(f['partition'])
        if len(labels) != node['components']:
            return self.fail(fid, 'partition size differs from the node')
        if kind == 'literature':
            _, name, g4 = f['source'].split(' ', 2)
            t = self.tables['g4'].get(name)
            if t is None:
                return self.fail(fid, f'{name} not in the tables')
            lo = int(t.strip('[]').split(';')[0])
            if t != g4 or lo < value:
                return self.fail(fid, f'table says {t}, fact says {g4} / >= {value}')
            if name == self.cert['target']:
                return self.fail(fid, "the target's own literature value")
            own, why_own = self.is_target_link(name)
            if own:
                return self.fail(fid, f"{name} is the target's own link ({why_own})")
            ok, why = prove_table_identity(self.node_gauss(f['node']), name, self.tables)
            if not ok:
                return self.fail(fid, f'node {f["node"]} not proved to be {name}: {why}')
            self.notes.append(f'fact {fid}: literature {name} {t} lower {lo}; identity: {why}')
            return True
        if kind == 'linking':
            if value != float('inf'):
                return self.fail(fid, 'a linking fact must say no surface exists')
            lk = self.node_gauss(f['node']).lk()
            blocks = defaultdict(list)
            for i, b in enumerate(labels):
                blocks[b].append(i)
            keys = sorted(blocks)
            for a in range(len(keys)):
                for b in range(a + 1, len(keys)):
                    if sum(lk[i][j] for i in blocks[keys[a]] for j in blocks[keys[b]]) != 0:
                        return True
            return self.fail(fid, 'the linking numbers allow this partition')
        if kind == 'witness':
            if f['from'] >= fid:
                return self.fail(fid, 'the fact read is not older')
            if not self.check_witness_edge(f):
                return False
            src = self.facts[f['from']]
            to_is_in = f['to_is_in']
            to, other = (f['in'], f['out']) if to_is_in else (f['out'], f['in'])
            if to != f['node'] or other != src['node']:
                return self.fail(fid, 'the fact and the one it read are not the edge\'s ends')
            to_map, other_map = (f['inMap'], f['outMap']) if to_is_in else (f['outMap'], f['inMap'])
            glued = [labels[to_map[j]] for j in range(len(to_map))]
            free, addition = glue(f['shape'], 'in' if to_is_in else 'out', glued, 0)
            other_labels = [None] * len(other_map)
            for i, v in enumerate(other_map):
                other_labels[v] = free[i]
            computed = normalize(other_labels)
            if computed != normalize(partition_labels(f['from_partition'])):
                return self.fail(fid, f'the cap gives partition {computed}, fact says '
                                      f'{f["from_partition"]}')
            stored = normalize(partition_labels(src['partition']))
            if not refines(computed, stored):
                return self.fail(fid, f'{computed} does not refine the partition read, {stored}')
            if addition != f['addition']:
                return self.fail(fid, f'the cap adds {addition}, fact says {f["addition"]}')
            src_value = self.value_of(src['value'])
            if src_value < self.value_of(f['from_value']):
                return self.fail(fid, 'the fact read is below the value used')
            if value > src_value - addition:
                return self.fail(fid, f'{value} > {src_value} - {addition}')
            return True
        return self.fail(fid, f'{kind}: replay not implemented')

    def run_lower(self):
        self.facts = {f['id']: f for f in self.cert['facts']}
        # The upper records a split fact subtracts, checked as in run().
        for rid in sorted(self.records):
            r = self.records[rid]
            try:
                ok = self.check_leaf(r) if r['kind'] == 'leaf' else \
                    self.check_witness(r) if r['kind'].startswith('witness') else self.check_split(r)
            except Exception as e:
                ok = self.fail(rid, f'{type(e).__name__}: {e}')
            print(f"  record {rid:>5} {r['kind']:<16} node {r['node']:>4} "
                  f"{r['partition']:<14} g{r['genus']}  {'ok' if ok else 'FAILED'}")
        for fid in sorted(self.facts):
            f = self.facts[fid]
            try:
                ok = self.check_lower_fact(f)
            except Exception as e:
                ok = self.fail(fid, f'{type(e).__name__}: {e}')
            print(f"  fact   {fid:>5} {f['kind']:<16} node {f['node']:>4} "
                  f"{f['partition']:<14} >= {f['value']}  {'ok' if ok else 'FAILED'}")
        top = self.facts[self.cert['top']]
        tnode = [n['id'] for n in self.cert['nodes'] if n['label'].startswith('target')][0]
        k = self.nodes[tnode]['components']
        goal = list(range(k)) if self.cert['goal'] == 'disjoint' else [0] * k
        if top['node'] != tnode or normalize(partition_labels(top['partition'])) != goal or \
                self.value_of(top['value']) < self.cert['goal_lower']:
            self.problems.append('the top fact is not the target meeting its lower goal')
        return not self.problems

    def run(self):
        ids = sorted(self.records)
        for rid in ids:
            r = self.records[rid]
            for c in r['children']:
                if c not in self.records or c >= rid:
                    self.fail(rid, f'child {c} missing or not older')
            try:
                if r['kind'] == 'leaf':
                    ok = self.check_leaf(r)
                elif r['kind'].startswith('witness'):
                    ok = self.check_witness(r)
                else:
                    ok = self.check_split(r)
            except Exception as e:
                ok = self.fail(rid, f'{type(e).__name__}: {e}')
            print(f"  record {rid:>5} {r['kind']:<16} node {r['node']:>4} "
                  f"{r['partition']:<14} g{r['genus']}  {'ok' if ok else 'FAILED'}")
        top = self.records[ids[-1]]
        goal_ok = top['genus'] <= self.cert['goal_genus'] and top['node'] == \
            [n['id'] for n in self.cert['nodes'] if n['label'].startswith('target')][0]
        if not goal_ok:
            self.problems.append('the last record is not the target meeting its goal')
        return not self.problems


def load_tables(knots, links, classes=None):
    """The tables, and (when the file exists) data/table_link_classes.csv:
    each table name to its class's canonical name. `classes` defaults to that
    file beside the link table."""
    t = {'g4': {}, 'pd': {}, 'classes': {}}
    for f in (knots, links):
        with open(f) as fh:
            for r in csv.DictReader(fh):
                t['g4'][r['Name']] = r['Genus-4D']
                t['pd'][r['Name']] = r[[k for k in r if k.startswith('PD')][0]]
    if classes is None:
        classes = os.path.join(os.path.dirname(links), 'table_link_classes.csv')
    if os.path.exists(classes):
        with open(classes) as fh:
            for r in csv.DictReader(fh):
                t['classes'][r['name']] = r['canonical']
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('certificate')
    ap.add_argument('--farsidediagram', default=DEFAULT_FSD)
    ap.add_argument('--knot-table', default=f'{ATLAS}/data/4d_smooth_slice_genus_13_crossings_pd_codes.csv')
    ap.add_argument('--link-table', default=f'{ATLAS}/data/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv')
    a = ap.parse_args()
    cert = json.load(open(a.certificate))
    ch = Checker(cert, a.farsidediagram, load_tables(a.knot_table, a.link_table))
    if 'facts' in cert:
        print(f"lower certificate: {cert['target']} genus >= {cert['lower']} "
              f"(goal {cert['goal_lower']}, {cert['goal']}), {len(cert['facts'])} facts, "
              f"{len(cert['records'])} records")
        ok = ch.run_lower()
    else:
        print(f"certificate: {cert['target']} genus <= {cert['genus']} "
              f"(goal {cert['goal_genus']}, {cert['goal']}), {len(cert['records'])} records")
        ok = ch.run()
    for n in ch.notes:
        print('  note:', n)
    for p in ch.problems:
        print('  PROBLEM:', p)
    print('VERDICT:', 'CERTIFIED' if ok else 'NOT CERTIFIED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()

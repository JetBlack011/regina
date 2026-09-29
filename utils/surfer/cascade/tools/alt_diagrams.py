#!/usr/bin/env python3
"""Alternative diagrams of a table knot, for a cascade target that its own
diagram leaves stuck: each is proved to be the same knot (up to mirror, which
g4 does not see) by an isometry of its exterior with the exterior of OUR
table PD (Gordon-Luecke), and written as a knotbuilder PD. The cascade's own
exact namer tags the target node with its table name as a second check.

A band a search cannot reach in one diagram's triangulation can be immediate
in another: 11a_164 resisted every hop shape on its table diagram (about
30,000 CPU-s) and fell in 10 s from the first alternative (2026-09-28).

usage: ~/.venvs/atlas/bin/python alt_diagrams.py <knot> <out.csv> <proofs.json>
           [--table CSV] [--max N]
out.csv: id,crossings,PD per line (ascending size, at most two per size)."""
import argparse, csv, json, random
import snappy  # noqa: F401 -- spherogram's exterior() needs it loaded
import spherogram

ap = argparse.ArgumentParser()
ap.add_argument('knot')
ap.add_argument('out')
ap.add_argument('proofs')
ap.add_argument('--table', default='/home/john/Projects/triangles/cobordism-atlas/data/'
                                   '4d_smooth_slice_genus_13_crossings_pd_codes.csv')
ap.add_argument('--max', type=int, default=8)
a = ap.parse_args()

row = next(r for r in csv.reader(open(a.table)) if r and r[0] == a.knot)
nums = [int(x) for x in row[1].replace('[', ' ').replace(']', ' ').replace(';', ' ').split()]
table = spherogram.Link([nums[i:i + 4] for i in range(0, len(nums), 4)])
ref = table.exterior()

random.seed(a.knot)
candidates = list(table.many_diagrams())
for _ in range(40):  # and some larger ones: random crossing-adding moves, lightly simplified
    L = table.copy()
    L.backtrack(random.randint(2, 6))
    L.simplify('basic')
    candidates.append(L)
seen, proved = set(), []
for L in candidates:
    code = L.PD_code()
    key = (len(code), str(L.DT_code()))
    if key in seen or len(L.link_components) != 1:
        continue
    seen.add(key)
    ext, iso = L.exterior(), False
    for _ in range(8):
        try:
            iso = bool(ext.is_isometric_to(ref))
            break
        except Exception:
            ext.randomize()
    if iso:
        proved.append({'crossings': len(code),
                       'pd': '[' + ';'.join('[' + ';'.join(map(str, c)) + ']' for c in code) + ']',
                       'proof': 'exterior isometric to the table PD exterior (SnapPy is_isometric_to)'})
proved.sort(key=lambda d: d['crossings'])
pick = []
for c in sorted({d['crossings'] for d in proved}):
    pick += [d for d in proved if d['crossings'] == c][:2]
pick = pick[:a.max]
with open(a.out, 'w', newline='') as fh:
    w = csv.writer(fh, lineterminator='\n')
    for i, d in enumerate(pick):
        w.writerow([f'd{i}', d['crossings'], d['pd']])
json.dump({'knot': a.knot, 'table_pd': row[1], 'diagrams': pick}, open(a.proofs, 'w'), indent=1)
print(f'{len(candidates)} candidates, {len(proved)} distinct and proved isometric; picked',
      [(f'd{i}', d['crossings']) for i, d in enumerate(pick)])

#!/usr/bin/env python3
"""Pool B: the rows the sweep did not settle (bounded, unresolved,
verified-assisted), each with a goal the cascade can prove constructively:
the literature value when the interval is a point, else its upper end.

Writes pool_B.csv in pool_A.csv's layout (name,g4,chain_rows,chain,
sweep_cpu_s,rows_with_provenance,pd), goal in g4, and prints a summary.
Read-only on the atlas."""
import csv, collections, re, sys

csv.field_size_limit(10**9)
A = '/home/john/Projects/triangles/cobordism-atlas'
out = sys.argv[1]

pd = {}
for f in ('data/4d_smooth_slice_genus_13_crossings_pd_codes.csv',
          'data/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv'):
    for r in csv.reader(open(f'{A}/{f}')):
        if r and r[0] != 'Name':
            pd[r[0]] = r[1]


def crossings(name):
    m = re.match(r'L?(\d+)', name.lstrip('m'))
    return int(m.group(1)) if m else 99


# Which searches a row had with a binary free of the row-drop bugs (D1/D2:
# the search side matched by name, orientation across surface components).
# Only 2026-09-c4's binary has the fixes; every earlier search may have
# discarded the row's surfaces, so its negative is weak.
fixed = {}
for p in csv.DictReader(open(f'{A}/results/search_provenance.csv')):
    if p['campaign'] == '2026-09-c4':
        fixed[p['row']] = f"c4:{p['new_witnesses']}w:{p['outcome']}"

rows = []
for r in csv.DictReader(open(f'{A}/results/verify_genus_v2.csv')):
    if r['status'] not in ('bounded', 'unresolved', 'verified-assisted'):
        continue
    lo, hi = r['literature_lo'], r['literature_hi']
    if hi in ('', 'inf', None) or r['knot'] not in pd:
        continue
    # A known answer the sweep did not reach: the literature is a point and
    # the atlas's constructive upper bound is missing or above it. (An open
    # interval whose upper end the sweep already has is research, not a
    # benchmark.)
    if lo != hi:
        continue
    goal = int(lo)
    if r['derived_hi'] not in ('', None) and int(r['derived_hi']) <= goal:
        continue
    rows.append(dict(name=r['knot'], status=r['status'], lo=lo, hi=hi, goal=goal,
                     dlo=r['derived_lo'], dhi=r['derived_hi'], x=crossings(r['knot']),
                     pd=pd[r['knot']], fixed=fixed.get(r['knot'], 'none')))
rows.sort(key=lambda r: (r['x'], r['status'], r['name']))
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, lineterminator='\n')
    w.writerow(['name', 'g4', 'chain_rows', 'chain', 'sweep_cpu_s', 'rows_with_provenance', 'pd',
                'status', 'literature', 'derived', 'fixed_search'])
    for r in rows:
        w.writerow([r['name'], r['goal'], 0, '', '', '', r['pd'], r['status'],
                    f"[{r['lo']};{r['hi']}]", f"[{r['dlo']};{r['dhi']}]", r['fixed']])
print(len(rows), 'rows;', collections.Counter(r['status'] for r in rows))
print('searched by a fixed binary (c4):', sum(r['fixed'] != 'none' for r in rows),
      '; only by binaries with the row-drop bugs:', sum(r['fixed'] == 'none' for r in rows))
for r in rows:
    if r['fixed'] != 'none':
        print('  fixed search:', r['name'], r['status'], r['fixed'])
print('by crossings:', sorted(collections.Counter(r['x'] for r in rows).items()))
print('knots', sum(not r['name'].startswith('L') for r in rows),
      'links', sum(r['name'].startswith('L') for r in rows))
print('point literature', sum(r['lo'] == r['hi'] for r in rows))
for r in rows[:12]:
    print(' ', r['name'], r['status'], f"lit [{r['lo']};{r['hi']}] derived [{r['dlo']};{r['dhi']}] goal {r['goal']}")

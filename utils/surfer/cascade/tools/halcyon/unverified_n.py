#!/usr/bin/env python3
"""The table knots of exactly N crossings the atlas has not verified --
including those with no status row at all -- as a pool CSV (pool_A.csv's
layout, goal = the literature value) for the cascade, with a summary.
Knots whose literature 4-genus is not a single value are listed and left out.

usage: unverified_n.py <out.csv> <N>"""
import collections, csv, re, sys

csv.field_size_limit(10**9)
A = '/home/john/Projects/triangles/cobordism-atlas'
out, N = sys.argv[1], int(sys.argv[2])
status = {r['knot']: r for r in csv.DictReader(open(f'{A}/results/verify_genus_v2.csv'))}
rows, counts, goals, odd = [], collections.Counter(), collections.Counter(), []
n = 0
for r in csv.reader(open(f'{A}/data/4d_smooth_slice_genus_13_crossings_pd_codes.csv')):
    if not r or r[0] == 'Name':
        continue
    m = re.fullmatch(r'(\d+)[an]?_(\d+)', r[0])
    if not m or int(m.group(1)) != N:
        continue
    n += 1
    s = status.get(r[0])
    st = s['status'] if s else 'no row'
    if st == 'verified':
        continue
    counts[st] += 1
    lit = r[2].strip()
    mm = re.fullmatch(r'\[?(\d+)(?:;(\d+))?\]?', lit)
    if not mm or (mm.group(2) is not None and mm.group(1) != mm.group(2)):
        odd.append((r[0], lit, st))
        continue
    goals[mm.group(1)] += 1
    rows.append([r[0], mm.group(1), 0, '', '', '', r[1], st, lit,
                 f"[{s['derived_lo']};{s['derived_hi']}]" if s else '', ''])
print(f'{n} table knots of {N} crossings; {sum(counts.values())} not verified:', dict(counts))
print('goals (literature g4):', dict(sorted(goals.items())))
print('left out, literature not a single value:', odd)
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, lineterminator='\n')
    w.writerow(['name', 'g4', 'chain_rows', 'chain', 'sweep_cpu_s', 'rows_with_provenance', 'pd',
                'status', 'literature', 'derived', 'fixed_search'])
    w.writerows(rows)
print(len(rows), 'targets written to', out)

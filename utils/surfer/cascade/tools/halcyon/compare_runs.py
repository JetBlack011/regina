#!/usr/bin/env python3
"""Pool A runs side by side: per target and in total, wall seconds (the
driver process, start to exit), hops and hop CPU, and whether the goal was met.

usage: compare_runs.py <tag>... (poolA = v1, poolA_v2, poolA_v3, ...)"""
import os, re, sys

R = os.path.expanduser('~/cascade-bench/results')
runs = {}
for tag in sys.argv[1:]:
    rows = {}
    path = os.path.join(R, tag + '.txt')
    if not os.path.exists(path):
        continue
    for line in open(path):
        m = re.match(r'(\S+) g4=(\d+) chain_rows=(\d+) .*hops=(\d+) cpu=(\d+) wall=(\d+) ', line)
        if m:
            rows[m.group(1)] = dict(g4=int(m.group(2)), chain=int(m.group(3)), hops=int(m.group(4)),
                                    cpu=int(m.group(5)), wall=int(m.group(6)),
                                    met='GOAL MET' in line, refused=line.count('refused'))
    runs[tag] = rows
tags = list(runs)
names = [n for n in runs[tags[0]]]
print(f'{"target":16s}' + ''.join(f'{t:>22s}' for t in tags))
for n in names:
    cells = []
    for t in tags:
        r = runs[t].get(n)
        cells.append(f'{r["wall"]:>5d}s {r["hops"]}h {r["cpu"]:>4d}cpu {"ok" if r["met"] else "--"}'
                     if r else ' ' * 22)
    print(f'{n:16s}' + ''.join(f'{c:>22s}' for c in cells))
for t in tags:
    rs = runs[t].values()
    print(f'{t}: {sum(r["met"] for r in rs)}/{len(rs)} met, {sum(r["hops"] for r in rs)} hops, '
          f'wall {sum(r["wall"] for r in rs)} s, hop CPU {sum(r["cpu"] for r in rs)} s, '
          f'refusals {sum(r["refused"] for r in rs)}')

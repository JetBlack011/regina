#!/usr/bin/env python3
"""The open pool: every table entry whose literature 4-genus is an interval
("[lo;hi]"), as cascade targets.

Writes one CSV row per entry:
  name,kind,crossings,components,lit_lo,lit_hi,goal_genus,master_witnesses,
  searched,pd
with goal_genus the interval's lower end (a proof of g4 <= lit_lo settles
the entry), master_witnesses the number of the master's witnesses with this
entry as subject, and searched whether any campaign's provenance row names
it. Ordered as the plan runs them: the [1;2] knots, then links by crossings
and components, then the [0;1] knots (the Dunfield-Gong sinks).

usage: open_pool.py <out.csv> [--atlas DIR] [--with-witnesses]

Read-only on the atlas. The master witness count streams cobordisms.csv's
subject column only, so it never holds the file in memory.
"""
import argparse
import collections
import csv
import os
import re
import sys

csv.field_size_limit(sys.maxsize)
KNOTS = 'data/4d_smooth_slice_genus_13_crossings_pd_codes.csv'
LINKS = 'data/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv'


def crossings(name):
    m = re.match(r'L?(\d+)', name)
    return int(m.group(1)) if m else 99


def components(name):
    m = re.search(r'\{([^}]*)\}', name)
    return 1 + len(m.group(1).split(';')) if m else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--atlas', default=os.environ.get(
        'CASCADE_ATLAS', '/home/john/Projects/triangles/cobordism-atlas'))
    ap.add_argument('--with-witnesses', action='store_true',
                    help='only entries the master holds witnesses for')
    a = ap.parse_args()

    entries = []
    for rel, kind in ((KNOTS, 'knot'), (LINKS, 'link')):
        with open(os.path.join(a.atlas, rel), newline='') as f:
            r = csv.reader(f)
            next(r)
            for row in r:
                m = re.match(r'^\[(\d+);(\d+)\]$', row[2].strip())
                if m:
                    entries.append({'name': row[0], 'kind': kind,
                                    'crossings': crossings(row[0]),
                                    'components': components(row[0]) if kind == 'link' else 1,
                                    'lit_lo': int(m.group(1)), 'lit_hi': int(m.group(2)),
                                    'pd': row[1]})
    names = {e['name'] for e in entries}

    witnesses = collections.Counter()
    with open(os.path.join(a.atlas, 'results/cobordisms.csv'), newline='') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) > 1 and row[1] in names:
                witnesses[row[1]] += 1
    searched = set()
    prov = os.path.join(a.atlas, 'results/search_provenance.csv')
    if os.path.exists(prov):
        with open(prov, newline='') as f:
            searched = {p['row'] for p in csv.DictReader(f)} & names

    def order(e):
        if e['kind'] == 'knot':
            return (0 if e['lit_lo'] >= 1 else 2, e['crossings'], e['name'])
        return (1, e['crossings'], e['components'], e['name'])

    rows = []
    for e in sorted(entries, key=order):
        e['goal_genus'] = e['lit_lo']
        e['master_witnesses'] = witnesses[e['name']]
        e['searched'] = int(e['name'] in searched)
        if a.with_witnesses and not e['master_witnesses']:
            continue
        rows.append(e)
    cols = ['name', 'kind', 'crossings', 'components', 'lit_lo', 'lit_hi',
            'goal_genus', 'master_witnesses', 'searched', 'pd']
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    kinds = collections.Counter((e['kind'], e['lit_lo'], e['lit_hi']) for e in rows)
    print(f'{len(rows)} open entries written to {a.out}: '
          + ', '.join(f'{k} [{lo};{hi}] {n}' for (k, lo, hi), n in sorted(kinds.items()))
          + f'; {sum(1 for e in rows if e["master_witnesses"])} with master witnesses, '
          f'{sum(e["searched"] for e in rows)} searched in a campaign')


if __name__ == '__main__':
    main()

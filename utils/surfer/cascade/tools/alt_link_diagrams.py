#!/usr/bin/env python3
"""Alternative diagrams of table LINKS, for cascade targets their own diagram
leaves stuck (alt_diagrams.py does knots, proving identity by exterior
isometry, which cannot pin a link: one complement belongs to many links).

Here identity holds by construction: every diagram comes from the table's
own PD by Reidemeister moves (spherogram's many_diagrams(), and random
backtracks lightly simplified), which keep the components and their
orientations. Each is checked to have the same number of components and the
same pairwise linking numbers, and at run time cascadesearch's exact namer
must name the target node as the entry itself, or it refuses the run.

usage: ~/.venvs/atlas/bin/python alt_link_diagrams.py <out.tsv> <name>... [--table CSV] [--max N]
out.tsv: name<TAB>label<TAB>pd per line (cascade_worker.py's names-file
form for an alternative diagram), at most two per crossing number, smallest
first, the table's own diagram excluded.
"""
import argparse
import random
import re

import snappy  # noqa: F401 -- spherogram's exterior() needs it loaded
import spherogram


def parse(pd):
    nums = [int(x) for x in re.findall(r'\d+', pd)]
    return [nums[i:i + 4] for i in range(0, len(nums), 4)]


def linking(L):
    m = L.linking_matrix()
    return sorted(m[i][j] for i in range(len(m)) for j in range(i + 1, len(m)))


def pd_text(L):
    return '[' + ';'.join('[' + ';'.join(map(str, c)) + ']' for c in L.PD_code()) + ']'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('names', nargs='+')
    ap.add_argument('--table', default='/home/john/Projects/triangles/cobordism-atlas/data/'
                                       'links_4d_smooth_slice_genus_11_crossings_pd_codes.csv')
    ap.add_argument('--max', type=int, default=6)
    a = ap.parse_args()
    table = {}
    with open(a.table) as f:
        next(f)
        for line in f:
            name, pd, _ = line.rstrip('\n').split(',', 2)
            table[name] = pd
    with open(a.out, 'w') as out:
        for name in a.names:
            L0 = spherogram.Link(parse(table[name]))
            n0, lk0 = len(L0.link_components), linking(L0)
            random.seed(name)
            cands = list(L0.many_diagrams())
            for _ in range(40):
                L = L0.copy()
                L.backtrack(random.randint(2, 6))
                L.simplify('basic')
                cands.append(L)
            seen = {pd_text(L0)}
            kept = []
            for L in cands:
                if len(L.link_components) != n0 or linking(L) != lk0:
                    continue  # never happens for Reidemeister moves; refused if it did
                text = pd_text(L)
                if text in seen:
                    continue
                seen.add(text)
                kept.append((len(L.PD_code()), text))
            kept.sort()
            pick = []
            for c in sorted({k for k, _ in kept}):
                pick += [t for k, t in kept if k == c][:2]
            pick = pick[:a.max]
            for i, text in enumerate(pick):
                out.write(f'{name}\td{i}\t{text}\n')
            print(f'{name}: {len(cands)} candidates, {len(kept)} distinct, wrote {len(pick)} '
                  f'({", ".join(str(len(parse(t)) // 1) for t in pick)} crossings)', flush=True)


if __name__ == '__main__':
    main()

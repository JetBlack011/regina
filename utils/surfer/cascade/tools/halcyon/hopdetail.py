#!/usr/bin/env python3
"""Per-hop detail of a Pool A run: where each hop's search wall goes beyond
its rounds and drain tail ("rest"), and which hops carry it.

usage: hopdetail.py <results dir, e.g. ~/cascade-bench/results/poolA_v5>"""
import json, os, statistics, sys

R = os.path.expanduser(sys.argv[1])
rows = []
for t in sorted(os.listdir(R)):
    p = os.path.join(R, t, 'cascade.jsonl')
    if not os.path.exists(p):
        continue
    for line in open(p):
        j = json.loads(line)
        if 'rounds' not in j:
            continue
        rounds = sum(j['rounds'])
        rest = j['search'] - rounds - j['drain_tail_s']
        extra = {k: j[k] for k in ('prototype', 'to_drain', 'seed_report') if k in j}
        rows.append(dict(rest=rest, drain=j['drain_tail_s'], search=j['search'], rounds=rounds,
                         t=t, hop=j['hop'], x=j['crossings'], w=j['witnesses'], **extra))
print(f"{len(rows)} hops: rest total {sum(r['rest'] for r in rows):.1f} s, "
      f"median {statistics.median(r['rest'] for r in rows):.2f} s; "
      f"drain tail median {statistics.median(r['drain'] for r in rows):.2f} s")
print(f"hops with rest > 1 s: {sum(r['rest'] > 1 for r in rows)}, "
      f"summing {sum(r['rest'] for r in rows if r['rest'] > 1):.1f} s")
for k in ('prototype', 'to_drain', 'seed_report'):
    if rows and k in rows[0]:
        print(f"{k}: total {sum(r[k] for r in rows):.1f} s, max {max(r[k] for r in rows):.2f} s")


def show(r):
    extra = ' '.join(f"{k} {r[k]:.2f}" for k in ('prototype', 'to_drain', 'seed_report') if k in r)
    return (f"  rest {r['rest']:.2f} drain {r['drain']:.2f} search {r['search']:.2f} "
            f"rounds {r['rounds']:.2f} {r['t']} hop{r['hop']} {r['x']}x {r['w']}w {extra}")


print('top 12 by rest:')
for r in sorted(rows, key=lambda r: -r['rest'])[:12]:
    print(show(r))
print('top 8 by drain tail:')
for r in sorted(rows, key=lambda r: -r['drain'])[:8]:
    print(show(r))

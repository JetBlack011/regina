#!/usr/bin/env python3
"""Per-hop phase timers from a Pool A run's cascade.jsonl files (v2 on):
the serial parts (row build and certification, search setup, adding kept
surfaces, naming new nodes, propagation) against the search itself.

usage: phases.py <results dir, e.g. ~/cascade-bench/results/poolA_v2> [target...]"""
import json, os, sys

R = os.path.expanduser(sys.argv[1])
targets = sys.argv[2:] or sorted(os.listdir(R))
tot = dict(row=0.0, setup=0.0, search=0.0, add=0.0, name_nodes=0.0, propagate=0.0, cpu=0.0, hops=0)
for t in targets:
    p = os.path.join(R, t, 'cascade.jsonl')
    if not os.path.exists(p):
        continue
    for line in open(p):
        j = json.loads(line)
        if 'row' not in j:
            continue
        if len(sys.argv) > 2:
            ser = j['row'] + j['setup'] + j['add'] + j['name_nodes'] + j['propagate']
            print(f"{t} hop {j['hop']} node {j['node']} ({j['crossings']}x, {j['witnesses']} kept):"
                  f" row {j['row']:.2f} setup {j['setup']:.2f} search {j['search']:.2f}"
                  f" add {j['add']:.2f} names {j['name_nodes']:.2f} prop {j['propagate']:.3f}"
                  f" | serial {ser:.2f} s, search {j['search']:.2f} s, search CPU {j['cpu']:.0f} s")
        for k in ('row', 'setup', 'search', 'add', 'name_nodes', 'propagate', 'cpu'):
            tot[k] += j[k]
        tot['hops'] += 1
serial = tot['row'] + tot['setup'] + tot['add'] + tot['name_nodes'] + tot['propagate']
wall = serial + tot['search']
# v4 on: the search's wall split into its IDDFS rounds and the drain tail
# (pending surfaces processed after the enumeration ended).
rounds, tail, tail_n, have = [], 0.0, 0, 0
for t in targets:
    p = os.path.join(R, t, 'cascade.jsonl')
    if not os.path.exists(p):
        continue
    for line in open(p):
        j = json.loads(line)
        if 'rounds' not in j:
            continue
        have += 1
        for i, r in enumerate(j['rounds']):
            if i == len(rounds):
                rounds.append(0.0)
            rounds[i] += r
        tail += j['drain_tail_s']
        tail_n += j['drain_tail']
if have:
    rest = tot['search'] - sum(rounds) - tail
    print(f"search wall {tot['search']:.1f} s over {have} hops = "
          + ' + '.join(f'round {i + 1} {r:.1f}' for i, r in enumerate(rounds))
          + f" + drain tail {tail:.1f} ({tail_n} surfaces) + rest {rest:.1f}")
if wall:
    print(f"{tot['hops']} hops: serial {serial:.1f} s ({serial / wall:.0%} of hop wall) = "
          f"row {tot['row']:.1f} + setup {tot['setup']:.1f} + add {tot['add']:.1f} + "
          f"names {tot['name_nodes']:.1f} + propagate {tot['propagate']:.1f}; "
          f"search {tot['search']:.1f} s wall, {tot['cpu']:.0f} s CPU "
          f"({tot['cpu'] / max(tot['search'], 1e-9):.1f} cores busy during it)")

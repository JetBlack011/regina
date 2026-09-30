#!/usr/bin/env python3
"""Where a cascade campaign's wall time goes, from each run's cascade.jsonl
(one line per hop, with phase timers) and its driver log.

Per hop: wall, cpu, setup, search (search + drain), drain_tail_s, naming_*_s,
add, name_nodes, propagate, assemble. Per run: the run's whole wall (from
.meta) versus the sum of its hops' walls (the rest is outside hops: startup,
master-witness loading, the store step, the lower report).

usage: phase_times.py ROOT [--since ISO] [--bin SUBSTR]"""
import collections, glob, json, os, re, sys

root = sys.argv[1]
since = sys.argv[sys.argv.index('--since') + 1] if '--since' in sys.argv else ''
binsub = sys.argv[sys.argv.index('--bin') + 1] if '--bin' in sys.argv else ''

tot = collections.Counter()
runs = 0
per_run = []
hop_rows = []
for meta in glob.glob(os.path.join(root, 'row_logs', '*.meta')):
    m = dict(l.rstrip('\n').split('=', 1) for l in open(meta) if '=' in l)
    if m.get('started_utc', '') < since:
        continue
    log = meta[:-5] + '.log'
    text = open(log).read() if os.path.exists(log) else ''
    if binsub and binsub not in text[:3000]:
        # the profile line names the binary only indirectly; check the meta
        pass
    work = meta[:-5] + '.cascade'
    hops = []
    jl = os.path.join(work, 'cascade.jsonl')
    if os.path.exists(jl):
        for l in open(jl):
            try:
                d = json.loads(l)
            except ValueError:
                continue
            if 'hop' in d and 'wall' in d:
                hops.append(d)
    if not hops:
        continue
    runs += 1
    wall = float(m.get('wall_s', 0) or 0)
    hw = sum(h.get('wall', 0) for h in hops)
    hc = sum(h.get('cpu', 0) for h in hops)
    store = re.search(r'witness store: .*signed in ([\d.]+) s', text)
    store_s = float(store[1]) if store else 0.0
    per_run.append((os.path.basename(meta)[:-5], wall, hw, hc, len(hops), store_s))
    tot['run_wall'] += wall
    tot['hop_wall'] += hw
    tot['hop_cpu'] += hc
    tot['store_sign'] += store_s
    for h in hops:
        for k in ('setup', 'search', 'drain_tail_s', 'naming_diagram_s', 'naming_fallback_s',
                  'naming_exact_s', 'add', 'name_nodes', 'propagate', 'assemble', 'row'):
            tot[k] += h.get(k, 0) or 0
        tot['hops'] += 1
        hop_rows.append(h)

print(f'{runs} runs, {tot["hops"]} hops')
print(f'run wall {tot["run_wall"]:.0f} s; inside hops {tot["hop_wall"]:.0f} s '
      f'({100 * tot["hop_wall"] / max(tot["run_wall"], 1):.0f}%); hop CPU {tot["hop_cpu"]:.0f} s '
      f'= {tot["hop_cpu"] / max(tot["hop_wall"], 1):.1f} cores busy on average inside hops')
print(f'store signing (after the hops) {tot["store_sign"]:.0f} s')
for k in ('row', 'setup', 'search', 'drain_tail_s', 'add', 'name_nodes', 'propagate', 'assemble'):
    print(f'  {k:14} {tot[k]:9.0f} s wall ({100 * tot[k] / max(tot["run_wall"], 1):5.1f}% of run wall)')
for k in ('naming_diagram_s', 'naming_fallback_s', 'naming_exact_s'):
    print(f'  {k:18} {tot[k]:9.0f} thread-s')
print('\nslowest runs (name, run wall, hop wall, hop cpu, hops, store sign):')
for r in sorted(per_run, key=lambda r: -r[1])[:8]:
    print('  %-28s %7.0f %7.0f %8.0f %4d %6.0f' % r)
print('\nslowest hops (wall, cpu, cores, search, drain, add, name_nodes, propagate, assemble, slowest name):')
for h in sorted(hop_rows, key=lambda h: -h.get('wall', 0))[:10]:
    cores = h.get('cpu', 0) / max(h.get('wall', 0), 1e-9)
    print('  %7.0f %8.0f %5.1f %7.0f %7.0f %6.0f %6.0f %6.0f %6.0f %6.1f' % (
        h.get('wall', 0), h.get('cpu', 0), cores, h.get('search', 0), h.get('drain_tail_s', 0),
        h.get('add', 0), h.get('name_nodes', 0), h.get('propagate', 0), h.get('assemble', 0),
        h.get('naming_slowest_s', 0)))

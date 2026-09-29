"""Summarise cascade runs: one row per run directory.

usage: python3 summarize_runs.py <runs dir> [--pool-a pool_A.csv]

Reads each run's cascade.jsonl (one line per hop) and driver.log. With
--pool-a, adds the sweep's recorded CPU for the same target (the rows along
its depends_on chain, from search_provenance.csv) and the ratio.
"""
import argparse, csv, glob, json, os


def run_row(d):
    hops = []
    p = os.path.join(d, 'cascade.jsonl')
    if os.path.exists(p):
        for line in open(p):
            try:
                hops.append(json.loads(line))
            except ValueError:
                pass
    log = open(os.path.join(d, 'driver.log')).read() if os.path.exists(
        os.path.join(d, 'driver.log')) else ''
    met = 'GOAL MET' in log
    basis = 'constructive' if '(constructive)' in log else (
        'assisted' if 'literature-assisted' in log else '')
    last = hops[-1] if hops else {}
    return {
        'run': os.path.basename(d),
        'met': met,
        'basis': basis,
        'hops': sum(1 for h in hops if 'cpu' in h),
        'refused': sum(1 for h in hops if 'refused' in h),
        'cpu_s': round(sum(h.get('cpu', 0) for h in hops)),
        'wall_s': round(sum(h.get('wall', 0) for h in hops)),
        'assemble_s': round(sum(h.get('assemble', 0) for h in hops)),
        'nodes': last.get('nodes', 0),
        'best': last.get('target_best'),
        'lower': last.get('target_lower'),
        'contradictions': last.get('contradictions', 0),
        'invariant_failures': last.get('invariant_failures', 0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs')
    ap.add_argument('--pool-a')
    a = ap.parse_args()
    sweep = {}
    if a.pool_a:
        for r in csv.DictReader(open(a.pool_a)):
            sweep[r['name']] = (int(r['sweep_cpu_s']), int(r['chain_rows']),
                                int(r['rows_with_provenance']), r['g4'])
    rows = []
    for d in sorted(glob.glob(os.path.join(a.runs, '*'))):
        if os.path.isdir(d) and os.path.exists(os.path.join(d, 'driver.log')):
            rows.append(run_row(d))
    cols = ['run', 'met', 'basis', 'hops', 'cpu_s', 'wall_s', 'assemble_s', 'nodes',
            'best', 'lower', 'contradictions', 'invariant_failures']
    if sweep:
        cols += ['g4', 'sweep_cpu_s', 'chain_rows', 'sweep/cascade']
        for r in rows:
            name = r['run'][2:] if r['run'].startswith('A_') else r['run']
            if name in sweep:
                s, n, k, g4 = sweep[name]
                r.update({'g4': g4, 'sweep_cpu_s': s, 'chain_rows': f'{k}/{n}',
                          'sweep/cascade': round(s / r['cpu_s'], 1) if r['cpu_s'] else ''})
    w = csv.DictWriter(open(os.devnull, 'w'), fieldnames=cols)
    print('\t'.join(cols))
    for r in rows:
        print('\t'.join(str(r.get(c, '')) for c in cols))
    met = [r for r in rows if r['met']]
    print(f'\n{len(met)}/{len(rows)} goals met; median CPU of successes '
          f'{sorted(r["cpu_s"] for r in met)[len(met) // 2] if met else "-"} s; '
          f'contradictions {sum(r["contradictions"] for r in rows)}, '
          f'invariant failures {sum(r["invariant_failures"] for r in rows)}')


if __name__ == '__main__':
    main()

"""Record CERTIFIED cascade proofs in the atlas's cascade store, with their
surfaces as ordinary witness rows.

usage: ~/.venvs/atlas/bin/python cascade_record.py --run <name> --store <dir>
           [--farsidediagram PATH] <target dir>...

A target dir is a cascadesearch --work directory holding certificate.json and
the check.txt that cascade_check.py wrote for it; only those whose check ends
"VERDICT: CERTIFIED" are recorded. The store (the atlas's results/cascade/)
gets, appended, never rewritten:

  proofs.csv     one row per proved target: goal, bound, hops, CPU, wall,
                 run, and where its certificate is
  witnesses.csv  every witness a proof uses, in cobordisms.csv's 13 columns,
                 its surface as a pair signature -- so the atlas's tools
                 (farsidediagram, farsidename, merge_cobordisms.py) read it as
                 they read the master
  proof_witnesses.csv  which proof and record each witness serves: a hop's
                 (in witnesses.csv), with the row it was found on and its
                 faces there, or the master's (--master-witnesses), by key
  nodes.csv      the links in those rows that are not table entries, named
                 cascade:<run>/<target>/n<id>, with their PD and Gauss data

How a surface becomes a row: its faces are rebuilt in the row's thickening
by farsidediagram --faces --pairsig, which first checks the thickening's
digest against the certificate's, then rebuilds the surface with the
search's own checks, and prints its pair signature over that thickening --
what verifyslicegenus would have recorded for it -- its genus, whether it is
connected, and its resolved vertices. The genus must be the certificate's,
and the signature must read back, through the route the atlas reads the
master by, as the same surface (genus, far-side curves, linking numbers).

Names: a node that is a table entry keeps its table name; a crossingless one
is `Unknot` or `<k>-component unlink`; any other is its cascade: name. A
split far side is its pieces joined by ` u `, unknots last, as the atlas
writes them.
"""
import argparse, csv, hashlib, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_FSD = os.path.normpath(os.path.join(HERE, '../../../../build/utils/surfer/farsidediagram'))
WITNESS_COLUMNS = ['kind', 'subject', 'subject_components', 'other', 'other_candidates',
                   'other_components', 'genus', 'tubed', 'pairsig', 'source_row', 'thicken_layers',
                   'max_faces', 'resolved_vertices']
PROOF_COLUMNS = ['run', 'target', 'goal', 'goal_genus', 'bound', 'hops', 'search_cpu_s', 'wall_s',
                 'witnesses', 'certificate']
LINK_COLUMNS = ['run', 'target', 'record', 'source', 'witness_key', 'subject', 'other', 'row_pd',
                'layers', 'build', 'faces']
NODE_COLUMNS = ['name', 'components', 'crossings', 'pd', 'signs', 'gauss', 'label']


def append(path, columns, rows):
    """Append rows (LF endings, as the C++ writes), with the header when new."""
    new = not os.path.exists(path)
    with open(path, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=columns, lineterminator='\n')
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)


def existing(path, key):
    if not os.path.exists(path):
        return set()
    with open(path, newline='') as fh:
        return {key(r) for r in csv.DictReader(fh)}


def node_name(run, target, node):
    if node.get('table'):
        return node['table']
    if node.get('label') == f'target {target}':
        return target  # a target beyond the tables keeps its own name (18nh_00000707)
    if not node.get('signs'):
        k = node['components']
        return 'Unknot' if k == 1 else f'{k}-component unlink'
    return f"cascade:{run}/{target}/n{node['id']}"


def far_name(names):
    unknots = sum(1 for n in names if n == 'Unknot')
    rest = sorted(n for n in names if n != 'Unknot')
    if not rest:
        return 'Unknot' if unknots == 1 else f'{unknots}-component unlink'
    return ' u '.join(rest + ['Unknot'] * unknots)


def sign(fsd, row_pd, layers, requests):
    """{record id: fields} for one row's surfaces, and the row's digest."""
    lines = ''.join(f"{rid} {','.join(str(f) for f in faces)}\n" for rid, faces in requests)
    out = subprocess.run([fsd, '--layers', str(layers), '--gauss', '--faces', '--pairsig', row_pd],
                         input=lines, capture_output=True, text=True, check=True).stdout.splitlines()
    row = dict(kv.split('=', 1) for kv in out[0].split(' ')[1:])
    got = {}
    for l in out[1:]:
        parts = l.split(' ')
        if parts[0] != 'W':
            continue
        if parts[2] != 'ok':
            raise RuntimeError(f'record {parts[1]}: farsidediagram: {l}')
        got[int(parts[1])] = dict(kv.split('=', 1) for kv in parts[3:])
    return row.get('build'), got


def read_back(fsd, row_pd, layers, sigs):
    """{id: fields} for pair signatures read as the atlas reads the master's
    (farsidediagram's pair-signature route: the surface carried onto the
    row's thickening by an isomorphism pinned by L x {0})."""
    lines = ''.join(f'{rid} {sig}\n' for rid, sig in sigs)
    out = subprocess.run([fsd, '--layers', str(layers), '--gauss', row_pd], input=lines,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    got = {}
    for l in out[1:]:
        parts = l.split(' ')
        if parts[0] == 'W':
            got[int(parts[1])] = (dict(kv.split('=', 1) for kv in parts[3:])
                                  if parts[2] == 'ok' else {'failed': l})
    return got


def lk_entries(text):
    return sorted(x for row in json.loads(text) for x in row)


def row_items(cert):
    """A proof's in-process witnesses grouped by the row they were found on:
    [((row PD, layers), [records])], in certificate order."""
    by_row = {}
    for r in cert['records']:
        if 'faces' in r:
            by_row.setdefault((r['row_pd'], r.get('layers', 2)), []).append(r)
    return list(by_row.items())


def sign_row(fsd, item):
    """One row's surfaces signed, and each signature read back."""
    (row_pd, layers), recs = item
    build, got = sign(fsd, row_pd, layers, [(r['id'], r['faces']) for r in recs])
    back = read_back(fsd, row_pd, layers, [(r['id'], got[r['id']]['pairsig']) for r in recs])
    return build, got, back


def record_target(d, run, max_faces, cert, items, results):
    """A proof's rows for the store, from its rows' signatures (`results`,
    one per item of row_items(cert))."""
    nodes = {n['id']: n for n in cert['nodes']}
    target = cert['target']
    witnesses, links, used = [], [], set()
    for ((row_pd, layers), recs), (build, got, back) in zip(items, results):
        for r in recs:
            f, b = got[r['id']], back.get(r['id'], {'failed': 'no line'})
            same = ('failed' not in b and b['genus'] == f['genus']
                    and b['components'] == f['components']
                    and lk_entries(b['lk']) == lk_entries(f['lk']))
            if not same:
                raise RuntimeError(f"{target} record {r['id']}: its pair signature does not read "
                                   f"back as the surface it was made from: {b}")
            if build != r['build']:
                raise RuntimeError(f"{target} record {r['id']}: the rebuilt thickening's digest "
                                   f"{build} is not the certificate's {r['build']}")
            f = got[r['id']]
            if 'shape' in r and int(f['genus']) != r['shape']['genus']:
                raise RuntimeError(f"{target} record {r['id']}: surface genus {f['genus']} is not "
                                   f"the certificate's {r['shape']['genus']}")
            subject_node = nodes[r['in']] if 'in' in r else nodes[r['node']]
            subject = node_name(run, target, subject_node)
            used.add(subject_node['id'])
            if r['kind'] == 'leaf':  # a direct witness: no far side
                kind, other, other_k = 'direct', '', 0
            else:
                names = [node_name(run, target, nodes[p['node']]) for p in r['pieces']]
                used.update(p['node'] for p in r['pieces'])
                kind, other = 'cobordism', far_name(names)
                other_k = sum(len(p['origins']) for p in r['pieces'])
            witnesses.append({
                'kind': kind, 'subject': subject,
                'subject_components': subject_node['components'], 'other': other,
                'other_candidates': other, 'other_components': other_k, 'genus': int(f['genus']),
                'tubed': 'false' if f['connected'] == '1' else 'true', 'pairsig': f['pairsig'],
                'source_row': subject, 'thicken_layers': layers, 'max_faces': max_faces,
                'resolved_vertices': '' if f['resolved'] == '0' else f['resolved']})
            links.append({'run': run, 'target': target, 'record': r['id'], 'source': 'hop',
                          'witness_key': hashlib.sha1(f['pairsig'].encode()).hexdigest()[:12],
                          'subject': subject, 'other': other, 'row_pd': row_pd, 'layers': layers,
                          'build': r['build'], 'faces': ','.join(str(x) for x in r['faces'])})
    # Witnesses the proof took from the master (--master-witnesses): already
    # atlas rows, so only referenced, by their key.
    for r in cert['records']:
        if 'pairsig' in r and 'faces' not in r:
            key = hashlib.sha1(r['pairsig'].encode()).hexdigest()[:12]
            if key != r.get('witness', '').split(':')[-1]:
                raise RuntimeError(f"{target} record {r['id']}: inline pair signature does not "
                                   f"hash to its key {r.get('witness')}")
            subject_node = nodes[r['in']]
            names = [node_name(run, target, nodes[p['node']]) for p in r.get('pieces', [])]
            used.add(subject_node['id'])
            used.update(p['node'] for p in r.get('pieces', []))
            links.append({'run': run, 'target': target, 'record': r['id'], 'source': 'master',
                          'witness_key': key, 'subject': node_name(run, target, subject_node),
                          'other': far_name(names) if names else '', 'row_pd': r.get('row_pd', ''),
                          'layers': r.get('layers', ''), 'build': '', 'faces': ''})
    node_rows = []
    for i in sorted(used):
        n = nodes[i]
        name = node_name(run, target, n)
        if name.startswith('cascade:'):
            node_rows.append({'name': name, 'components': n['components'],
                              'crossings': len(n.get('signs', [])), 'pd': n.get('pd', ''),
                              'signs': json.dumps(n.get('signs', [])),
                              'gauss': json.dumps(n.get('gauss', [])), 'label': n.get('label', '')})
    hops = [json.loads(l) for l in open(os.path.join(d, 'cascade.jsonl')) if l.strip()]
    proof = {'run': run, 'target': target, 'goal': cert.get('goal', ''),
             'goal_genus': cert.get('goal_genus', ''), 'bound': cert.get('genus', ''),
             'hops': len(hops), 'search_cpu_s': round(sum(h.get('cpu', 0) for h in hops)),
             'wall_s': round(sum(h.get('wall', 0) for h in hops)), 'witnesses': len(witnesses),
             'certificate': f'{run}/{os.path.basename(os.path.normpath(d))}/certificate.json'}
    return proof, witnesses, links, node_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--store', required=True)
    ap.add_argument('--farsidediagram', default=DEFAULT_FSD)
    ap.add_argument('--max-faces', type=int, default=5)
    ap.add_argument('--jobs', type=int, default=4, help='rows signed at once')
    ap.add_argument('dirs', nargs='+')
    a = ap.parse_args()
    os.makedirs(a.store, exist_ok=True)
    csv.field_size_limit(10**9)
    done = existing(os.path.join(a.store, 'proofs.csv'), lambda r: (r['run'], r['target']))
    known_nodes = existing(os.path.join(a.store, 'nodes.csv'), lambda r: r['name'])
    # One row per surface: a proof may use one witness twice (once in each
    # direction), and proofs may share witnesses.
    known_sigs = existing(os.path.join(a.store, 'witnesses.csv'), lambda r: r['pairsig'])
    todo = []
    for d in a.dirs:
        check = os.path.join(d, 'check.txt')
        if not os.path.exists(check) or not open(check).read().rstrip().endswith('VERDICT: CERTIFIED'):
            print(f'skipped (not CERTIFIED): {d}')
            continue
        cert = json.load(open(os.path.join(d, 'certificate.json')))
        if (a.run, cert['target']) in done:
            print(f"already recorded: {a.run} {cert['target']}")
            continue
        todo.append((d, cert, row_items(cert)))
    # Every row of every proof signed on one pool: most proofs have one or two
    # rows, so signing proof by proof leaves the cores idle.
    tasks = [(k, i) for k, (_, _, items) in enumerate(todo) for i in range(len(items))]
    results = {}
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        futures = {pool.submit(sign_row, a.farsidediagram, todo[k][2][i]): (k, i) for k, i in tasks}
        for fut in futures:
            k, i = futures[fut]
            try:
                results[(k, i)] = fut.result()
            except Exception as e:
                results[(k, i)] = e
    for k, (d, cert, items) in enumerate(todo):
        failed = [results[(k, i)] for i in range(len(items)) if isinstance(results[(k, i)], Exception)]
        try:
            if failed:
                raise failed[0]
            proof, witnesses, links, node_rows = record_target(
                d, a.run, a.max_faces, cert, items, [results[(k, i)] for i in range(len(items))])
        except Exception as e:
            print(f"FAILED {cert['target']}: {e}")
            continue
        target = cert['target']
        fresh = []
        for w in witnesses:
            if w['pairsig'] not in known_sigs:
                known_sigs.add(w['pairsig'])
                fresh.append(w)
        append(os.path.join(a.store, 'witnesses.csv'), WITNESS_COLUMNS, fresh)
        append(os.path.join(a.store, 'proof_witnesses.csv'), LINK_COLUMNS, links)
        append(os.path.join(a.store, 'nodes.csv'), NODE_COLUMNS,
               [n for n in node_rows if n['name'] not in known_nodes])
        known_nodes.update(n['name'] for n in node_rows)
        append(os.path.join(a.store, 'proofs.csv'), PROOF_COLUMNS, [proof])
        print(f"recorded {target}: bound {proof['bound']}, {len(witnesses)} witnesses "
              f"({len(fresh)} new), {len(node_rows)} cascade nodes")


if __name__ == '__main__':
    sys.exit(main())

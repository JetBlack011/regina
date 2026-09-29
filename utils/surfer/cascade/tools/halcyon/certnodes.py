#!/usr/bin/env python3
"""The nodes a certificate's proof records use: id, components, crossings,
table name (when the cascade named one), and each record's source.

usage: certnodes.py <certificate.json>"""
import json, sys

c = json.load(open(sys.argv[1]))
nodes = {n['id']: n for n in c['nodes']}
used = []
for r in c['records']:
    used.append(r['node'])
    extra = {k: r[k] for k in ('source', 'kind', 'genus', 'partition') if k in r}
    print('record', r['id'], extra)
print('nodes used:')
for i in dict.fromkeys(used):
    n = nodes[i]
    print(' ', i, {k: n[k] for k in n if k not in ('gauss', 'signs', 'pd')},
          'crossings', len(n.get('signs', [])))

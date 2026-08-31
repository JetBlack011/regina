#!/usr/bin/env python3
"""Turn verifyslicegenus's witness file into a self-contained, double-clickable
visualization of the cobordism graph:

    python3 export_cobordism_graph.py

Companion to export_graph_data.py, which does the same job for a single
surfer run over one knot. This one covers the whole accumulated graph that
verifyslicegenus builds across runs -- every cobordism witness ever recorded,
and what each one grounds.

Pure data-shaping over cobordisms.csv / verify_genus_v2.csv; no regina, no
snappy, no network access. Layout is precomputed here and the result spliced
into cobordism_graph_template.html.

WHY RADIAL BY DISTANCE FROM THE UNKNOT. With --no-cone there are no direct
witnesses, so the only constructive anchors are the unknot (bounds a disc) and
the n-component unlinks (bound n disjoint discs). Every verified genus is
therefore the end of a chain of cobordisms back to one of those. Distance from
the Unknot is thus not an arbitrary aesthetic choice -- it is the quantity the
whole method turns on, and the picture it gives (a verified core, an
unverified rim, and a detached cloud that no chain reaches) is the actual
state of the computation.
"""
import argparse
import collections
import csv
import json
import math
import os
import re
import sys

csv.field_size_limit(sys.maxsize)

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA = os.path.normpath(
    os.path.join(HERE, "..", "..", "..", "build", "utils", "surfer"))
TEMPLATE = os.path.join(HERE, "cobordism_graph_template.html")
PLACEHOLDER = "__COBORDISM_DATA_JSON__"


def kind_of(name):
    """Node class. Drives colour and shape, and matters mathematically: an
    unlink is an ANCHOR (chi* = n, no component penalty), while an unnamed
    multi-component node is inert -- a complement we could not name, whose
    component orientations are therefore unknowable."""
    if name == "Unknot":
        return "unknot"
    if name.endswith("-component unlink"):
        return "unlink"
    if re.match(r"^L\d+[an]", name):
        return "link"
    # Superscript notation: the superscript IS the component count, so
    # 10^2_102 is a 2-component LINK, not a knot. Reading these as knots
    # inflated the knot count from 165 to 246.
    if re.match(r"^\d+\^\d+_", name):
        return "link"
    if re.match(r"^\d+_\d+$", name):
        return "knot"
    return "unnamed"


def components_of(name):
    m = re.match(r"^(\d+)-component unlink$", name)
    if m:
        return int(m.group(1))
    m = re.match(r"^\d+\^(\d+)_", name)          # 10^2_102 -> 2 components
    if m:
        return int(m.group(1))
    m = re.search(r"\{(.*)\}$", name)
    return len(m.group(1).split(";")) + 1 if m else 1


def load(data_dir):
    cob = os.path.join(data_dir, "cobordisms.csv")
    out = os.path.join(data_dir, "verify_genus_v2.csv")
    with open(cob, newline="") as f:
        edges = list(csv.DictReader(f))
    rows = {}
    if os.path.exists(out):
        with open(out, newline="") as f:
            rows = {r["knot"]: r for r in csv.DictReader(f)}
    return edges, rows


def build(edges, rows):
    adj = collections.defaultdict(set)
    for e in edges:
        a, b = e["subject"], e["other"]
        adj[a].add(b)
        adj[b].add(a)
    names = sorted(adj)

    # BFS from the Unknot, recording a parent so the layout can place children
    # near their parent and the UI can show the chain that grounds a node.
    dist, parent = {}, {}
    if "Unknot" in adj:
        dist["Unknot"] = 0
        queue = collections.deque(["Unknot"])
        while queue:
            cur = queue.popleft()
            for nb in sorted(adj[cur]):
                if nb not in dist:
                    dist[nb] = dist[cur] + 1
                    parent[nb] = cur
                    queue.append(nb)
    unreachable = [n for n in names if n not in dist]
    max_d = max(dist.values()) if dist else 0
    detached_ring = max_d + 2          # visibly separated from the reachable rings
    for n in unreachable:
        dist[n] = detached_ring

    # Angles: within each ring, order by the parent's angle so subtrees stay
    # contiguous; ties broken by status then name so colour forms arcs rather
    # than speckle.
    STATUS_ORDER = {"verified": 0, "verified-assisted": 1, "improved": 2,
                    "pinned": 3, "bounded": 4, "unresolved": 5}

    def status_of(n):
        return rows.get(n, {}).get("status", "untracked")

    rings = collections.defaultdict(list)
    for n in names:
        rings[dist[n]].append(n)

    angle = {}
    for d in sorted(rings):
        members = rings[d]
        if d == 0:
            angle[members[0]] = 0.0
            continue
        members.sort(key=lambda n: (
            angle.get(parent.get(n), 0.0),
            STATUS_ORDER.get(status_of(n), 9),
            n))
        for i, n in enumerate(members):
            angle[n] = 2 * math.pi * i / len(members)

    # Ring radii, with alternate nodes staggered inward/outward so the dense
    # rings (215 nodes at distance 4) do not collide.
    R0, DR, STAGGER = 0.0, 150.0, 22.0
    coords = {"Unknot": (0.0, 0.0)} if "Unknot" in adj else {}
    for d, members in rings.items():
        if d == 0:
            continue
        for i, n in enumerate(members):
            r = R0 + DR * d + (STAGGER if i % 2 else -STAGGER)
            coords[n] = (r * math.cos(angle[n]), r * math.sin(angle[n]))

    index = {n: i for i, n in enumerate(names)}
    node_json = []
    for n in names:
        r = rows.get(n, {})
        x, y = coords[n]
        node_json.append({
            "n": n,
            "x": round(x, 1),
            "y": round(y, 1),
            "k": kind_of(n),
            "s": status_of(n),
            "d": dist[n] if n not in unreachable else -1,
            "c": components_of(n),
            "lo": r.get("literature_lo", ""),
            "hi": r.get("literature_hi", ""),
            "dh": r.get("derived_hi", ""),
            "b": r.get("witness_basis", ""),
            "via": r.get("via_knot", ""),
            "p": index[parent[n]] if n in parent else -1,
            "deg": len(adj[n]),
        })

    edge_json = []
    seen = set()
    for e in edges:
        a, b = e["subject"], e["other"]
        g, tb = int(e["genus"]), e["tubed"] == "true"
        key = (a, b, g, tb)
        if key in seen:
            continue
        seen.add(key)
        edge_json.append({"s": index[a], "t": index[b], "g": g,
                          "tb": 1 if tb else 0})

    counts = collections.Counter(n["s"] for n in node_json)
    kinds = collections.Counter(n["k"] for n in node_json)
    return {
        "nodes": node_json,
        "edges": edge_json,
        "meta": {
            "nodes": len(node_json),
            "edges": len(edge_json),
            "witnesses": len(edges),
            "maxRing": detached_ring,
            "status": dict(counts),
            "kinds": dict(kinds),
            "tubed": sum(1 for e in edge_json if e["tb"]),
            "unreachable": len(unreachable),
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=DEFAULT_DATA,
                    help="directory holding cobordisms.csv and verify_genus_v2.csv")
    ap.add_argument("--out", default=os.path.join(HERE, "cobordism_atlas.html"))
    a = ap.parse_args()

    edges, rows = load(a.data_dir)
    data = build(edges, rows)

    with open(TEMPLATE) as f:
        html = f.read()
    if PLACEHOLDER not in html:
        sys.exit(f"template {TEMPLATE} is missing {PLACEHOLDER}")
    html = html.replace(PLACEHOLDER, json.dumps(data, separators=(",", ":")))
    with open(a.out, "w") as f:
        f.write(html)

    m = data["meta"]
    print(f"{m['nodes']} nodes, {m['edges']} distinct edges "
          f"({m['witnesses']} witnesses, {m['tubed']} tubed)")
    print(f"  status: {m['status']}")
    print(f"  kinds:  {m['kinds']}")
    print(f"  unreachable from the Unknot: {m['unreachable']}")
    print(f"wrote {a.out} ({os.path.getsize(a.out)//1024} KB)")


if __name__ == "__main__":
    main()

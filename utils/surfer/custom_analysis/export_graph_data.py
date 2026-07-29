#!/usr/bin/env python
"""Turn concordances.json into a compact, pre-laid-out graph, then splice it
into graph_viz_template.html to produce the self-contained, double-click-
able concordance_atlas.html -- the whole visualization build in one command:

    python3 export_graph_data.py

Run under plain python3 (no snappy/regina needed -- this is pure
data-shaping over already-computed concordances.json).
"""
import argparse
import json
import math
from collections import defaultdict

CONCORDANCES_JSON = "concordances.json"
OUTPUT_JSON = "graph_data.json"
TEMPLATE_HTML = "graph_viz_template.html"
OUTPUT_HTML = "concordance_atlas.html"
PLACEHOLDER = "__GRAPH_DATA_JSON__"
HUB_NAME = "K13n65"

# Surface types found in the data, ordered by frequency (most common first) --
# this fixes the categorical color-slot assignment used by the artifact.
SURFACE_TYPE_ORDER = [
    (True, 0, 3),
    (False, 1, 3),
    (True, 0, 4),
    (False, 1, 2),
    (False, 2, 2),
    (True, 1, 2),
    (True, 0, 2),
]


def surface_key(s):
    return (s["orientable"], s["genus"], s["punctures"])


def slice_status(genus):
    if genus is None:
        return "UNKNOWN"
    if isinstance(genus, list):
        lo, hi = genus
        return "OPEN" if lo <= 0 <= hi else "CONFIRMED_NOT_SLICE"
    if genus == 0:
        return "CONFIRMED_SLICE"
    return "CONFIRMED_NOT_SLICE"


def build_graph(concordances):
    nodes = {}  # name -> {genus, slice_status}
    edges = []  # (a, b, [surface_type_index,...])
    adjacency = defaultdict(set)

    surface_index = {s: i for i, s in enumerate(SURFACE_TYPE_ORDER)}

    for entry in concordances:
        a = entry["links"][0]["name"]
        b = entry["links"][1]["name"]
        ga = entry["links"][0]["4d_genus"]
        gb = entry["links"][1]["4d_genus"]
        for name, genus in ((a, ga), (b, gb)):
            if name not in nodes:
                nodes[name] = {"genus": genus, "slice_status": slice_status(genus)}
        types = sorted({surface_index[surface_key(s)] for s in entry["concordances"]})
        edges.append({"source": a, "target": b, "types": types})
        adjacency[a].add(b)
        if b != a:
            adjacency[b].add(a)

    degree = {name: len(neigh) for name, neigh in adjacency.items()}
    for name in nodes:
        degree.setdefault(name, 0)

    return nodes, edges, degree, adjacency


def dominant_surface_type(name, edges_by_node):
    """For a leaf node, the surface type of its one edge (for angular bucketing)."""
    es = edges_by_node.get(name, [])
    if not es:
        return 0
    return min(es[0]["types"]) if es[0]["types"] else 0


def layout(nodes, degree, edges, hub_ring_count=40):
    edges_by_node = defaultdict(list)
    for e in edges:
        edges_by_node[e["source"]].append(e)
        if e["target"] != e["source"]:
            edges_by_node[e["target"]].append(e)

    others = [n for n in nodes if n != HUB_NAME]
    others.sort(key=lambda n: -degree.get(n, 0))

    hub_ring_names = set(others[:hub_ring_count])
    leaves = [n for n in others if n not in hub_ring_names]

    coords = {HUB_NAME: (0.0, 0.0)}

    # Middle ring: secondary hubs, grouped by slice_status then sorted by degree.
    status_order = {"CONFIRMED_SLICE": 0, "OPEN": 1, "CONFIRMED_NOT_SLICE": 2, "UNKNOWN": 3}
    hub_ring = sorted(
        hub_ring_names,
        key=lambda n: (status_order.get(nodes[n]["slice_status"], 9), -degree.get(n, 0)),
    )
    # Staggered across two sub-radii (alternating by index) so adjacent labels
    # in angle don't also sit at the same radius -- doubles the effective
    # spacing between the 40 labeled nodes without shrinking the font.
    R_MID, R_MID_STAGGER = 400.0, 55.0
    n_mid = max(len(hub_ring), 1)
    for i, name in enumerate(hub_ring):
        angle = 2 * math.pi * i / n_mid
        r = R_MID + (R_MID_STAGGER if i % 2 else -R_MID_STAGGER)
        coords[name] = (r * math.cos(angle), r * math.sin(angle))

    # Outer ring: leaves, angularly bucketed by (slice_status, dominant surface type)
    # so same-colored/same-edge-style nodes cluster into contiguous arcs.
    leaves.sort(
        key=lambda n: (
            status_order.get(nodes[n]["slice_status"], 9),
            dominant_surface_type(n, edges_by_node),
            n,
        )
    )
    R_OUT = 800.0
    n_out = max(len(leaves), 1)
    for i, name in enumerate(leaves):
        angle = 2 * math.pi * i / n_out
        # small radius jitter by degree so degree>1 leaves (rare) sit slightly
        # inward, degree==1 pure leaves sit on the outer rim
        r = R_OUT - (10.0 if degree.get(name, 0) > 1 else 0.0)
        coords[name] = (r * math.cos(angle), r * math.sin(angle))

    # The isolated component (not connected to the hub at all): place off to
    # the side rather than at the origin-adjacent coordinate space.
    placed = set(coords)
    stray = [n for n in nodes if n not in placed]
    for i, name in enumerate(stray):
        coords[name] = (R_OUT + 120.0, i * 40.0 - (len(stray) - 1) * 20.0)

    return coords, hub_ring_names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=CONCORDANCES_JSON)
    ap.add_argument("--output", default=OUTPUT_JSON, help="intermediate graph data JSON")
    ap.add_argument("--template", default=TEMPLATE_HTML)
    ap.add_argument("--html-output", default=OUTPUT_HTML)
    ap.add_argument("--data-only", action="store_true", help="skip the HTML build, just write the JSON")
    args = ap.parse_args()

    with open(args.input) as f:
        concordances = json.load(f)

    nodes, edges, degree, adjacency = build_graph(concordances)
    coords, hub_ring_names = layout(nodes, degree, edges)

    out_nodes = []
    for name, info in nodes.items():
        x, y = coords[name]
        out_nodes.append({
            "id": name,
            "x": round(x, 2),
            "y": round(y, 2),
            "genus": info["genus"],
            "status": info["slice_status"],
            "degree": degree.get(name, 0),
            "hub": name == HUB_NAME,
            "secondary_hub": name in hub_ring_names,
        })

    out = {
        "hub": HUB_NAME,
        "surface_types": [
            {"orientable": o, "genus": g, "punctures": p} for (o, g, p) in SURFACE_TYPE_ORDER
        ],
        "nodes": out_nodes,
        "edges": edges,
    }

    data_json = json.dumps(out, separators=(",", ":"))
    with open(args.output, "w") as f:
        f.write(data_json)

    status_counts = defaultdict(int)
    for n in out_nodes:
        status_counts[n["status"]] += 1
    print(f"nodes: {len(out_nodes)}  edges: {len(edges)}")
    print(f"status counts: {dict(status_counts)}")
    print(f"secondary hubs: {len(hub_ring_names)}")
    print(f"wrote {args.output}")

    if args.data_only:
        return

    with open(args.template) as f:
        template = f.read()
    if PLACEHOLDER not in template:
        raise ValueError(f"{args.template} is missing the {PLACEHOLDER} placeholder")
    with open(args.html_output, "w") as f:
        f.write(template.replace(PLACEHOLDER, data_json))
    print(f"wrote {args.html_output}")


if __name__ == "__main__":
    main()

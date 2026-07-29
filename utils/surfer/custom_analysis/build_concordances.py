#!/usr/bin/env python
"""Collapse the row-per-surface identified CSV down to a list of
concordances: distinct (knot/link, knot/link) pairs that are cobounded by
some surface in the dataset, along with every distinct surface type that
connects them and (where known) each side's smooth 4D slice genus.

Input: the output of identify_boundaries.py
       (13n_65_2_layer_collar_identified.csv), plus
       4d_smooth_slice_genus_13_crossings.csv.
Output: concordances.json

Run under .venv/bin/python. Pure CSV/JSON post-processing -- no snappy or
regina needed here.

See /home/john/.claude/plans/you-ll-notice-that-you-re-zany-meadow.md for
the design rationale.
"""
import argparse
import csv
import json
import re
import sys

csv.field_size_limit(sys.maxsize)

INPUT_CSV = "13n_65_2_layer_collar_identified.csv"
GENUS_CSV = "4d_smooth_slice_genus_13_crossings.csv"
OUTPUT_JSON = "concordances.json"

_KNOTINFO_RE = re.compile(r"^K(\d+)([an])(\d+)$")


def load_genus_table(path):
    table = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["Name"]
            raw = row["Genus-4D"]
            m = re.match(r"^\[(\d+);(\d+)\]$", raw)
            table[name] = [int(m.group(1)), int(m.group(2))] if m else int(raw)
    return table


def lookup_4d_genus(alias_names, genus_table):
    for name in alias_names:
        if name in genus_table:
            return genus_table[name]
        m = _KNOTINFO_RE.match(name)
        if m:
            normalized = f"{m.group(1)}{m.group(2)}_{m.group(3)}"
            if normalized in genus_table:
                return genus_table[normalized]
    return None


def component_identity(entry):
    """Return (name, alias_names, is_knot) for one boundary component's
    detail entry (one value of the per-row boundary_identification_json
    dict). A component is a "knot" only when it's a single piece -- a
    `whole` link (2+ pieces) is never a knot, regardless of whether it
    resolved."""
    if entry["whole"] is not None:
        w = entry["whole"]
        if w["kind"] == "identified":
            return w["best_name"], w["names"], False
        return w["token"], [], False
    piece = entry["pieces"][0]
    if piece["kind"] == "unknot":
        return "Unknot", [], True
    if piece["kind"] == "identified":
        return piece["best_name"], piece["names"], True
    return piece["token"], [], True


def process(input_path, genus_table, limit=None):
    pairs = {}  # (name_a, name_b) -> set of (orientable, genus, punctures)
    identities = {}  # name -> (alias_names, is_knot), first-seen wins
    n = 0
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n += 1
            detail = json.loads(row["boundary_identification_json"])
            name1, aliases1, isk1 = component_identity(detail["1"])
            name2, aliases2, isk2 = component_identity(detail["2"])
            identities.setdefault(name1, (aliases1, isk1))
            identities.setdefault(name2, (aliases2, isk2))

            surface = (
                row["orientable"].strip().lower() == "true",
                int(row["genus"]),
                int(row["punctures"]),
            )
            key = tuple(sorted((name1, name2)))
            pairs.setdefault(key, set()).add(surface)

            if limit and n >= limit:
                break
    return pairs, identities, n


def build_output(pairs, identities, genus_table):
    genus_cache = {}
    fallback_count = 0
    knot_miss_count = 0

    def genus_for(name, aliases):
        if name not in genus_cache:
            if not aliases:
                genus_cache[name] = None
            else:
                genus_cache[name] = lookup_4d_genus(aliases, genus_table)
        return genus_cache[name]

    out = []
    for (name_a, name_b), surfaces in pairs.items():
        aliases_a, isk_a = identities[name_a]
        aliases_b, isk_b = identities[name_b]
        g_a = genus_for(name_a, aliases_a) if name_a != "Unknot" else 0
        g_b = genus_for(name_b, aliases_b) if name_b != "Unknot" else 0
        for name, aliases, is_knot, g in (
            (name_a, aliases_a, isk_a, g_a), (name_b, aliases_b, isk_b, g_b),
        ):
            if is_knot and not aliases and name != "Unknot":
                fallback_count += 1
            elif is_knot and aliases and g is None:
                knot_miss_count += 1
        concordances = [
            {"orientable": o, "genus": gg, "punctures": p}
            for o, gg, p in sorted(surfaces)
        ]
        out.append({
            "links": [
                {"name": name_a, "4d_genus": g_a},
                {"name": name_b, "4d_genus": g_b},
            ],
            "concordances": concordances,
        })
    return out, fallback_count, knot_miss_count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=INPUT_CSV)
    ap.add_argument("--genus-csv", default=GENUS_CSV)
    ap.add_argument("--output", default=OUTPUT_JSON)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    genus_table = load_genus_table(args.genus_csv)
    print(f"loaded {len(genus_table)} genus entries", file=sys.stderr)

    pairs, identities, n_rows = process(args.input, genus_table, limit=args.limit)
    print(f"processed {n_rows} rows -> {len(pairs)} unique pairs, {len(identities)} unique identities", file=sys.stderr)

    out, fallback_count, knot_miss_count = build_output(pairs, identities, genus_table)
    print(f"identities using raw-isosig fallback (unresolved): {fallback_count}", file=sys.stderr)
    print(f"identified knot names with no genus-table hit: {knot_miss_count}", file=sys.stderr)

    with open(args.output, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {len(out)} concordance entries -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

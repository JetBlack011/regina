#!/usr/bin/env python3
"""
Validates ../knot_census_names.csv (produced by gen_knot_census_names.cpp)
against SnapPy, then emits ../rolfsentable.h: a census-name -> Rolfsen-name
lookup table for linkcomplement.cpp.

Usage (run with the venv that has SnapPy installed):
    utils/surfer/.venv/bin/python3 gen_knot_names_header.py

The table is keyed by *base* manifold name (e.g. "m004"), not the full raw
hit string Census::lookup() returns (e.g. "m004 : #1"). The " : #N" suffix
disambiguates which of possibly several combinatorially distinct ideal
triangulations of the *same* manifold was matched -- Regina's
cusped-hyp-or-census-9 stores one entry per triangulation discovered during
census construction, and many manifolds have several (confirmed: every
group of entries sharing a base name has identical hyperbolic volume).
Keying on the full string would only ever recognize whichever one specific
triangulation knotbuilder.cpp + simplify() happened to produce for a given
PD code, and silently miss the same manifold reached via any other
triangulation (e.g. one a live surfer search stumbles on mid-simplification).

Validation performed:
  - Every row hitting "Cusped hyperbolic census (...)" is cross-checked
    against snappy.Manifold(rolfsen_name).identify(): the bare census name
    (e.g. "m004" from "m004 : #1") must appear among SnapPy's independent
    identification of the same knot.
  - Every Rolfsen knot with *no* census hit is cross-checked as genuinely
    non-hyperbolic (near-zero volume) via SnapPy. A "no hit" knot that
    SnapPy says *is* hyperbolic is an anomaly -- it just means the
    complement needs more than 9 tetrahedra to fall inside Regina's shipped
    census, not a bug, but it's worth seeing called out explicitly.
  - Collisions (the same base census name mapping to two different Rolfsen
    names) are detected and excluded from the generated table rather than
    silently resolved.
"""

import csv
import pathlib
import sys

import snappy

HERE = pathlib.Path(__file__).resolve().parent
CSV_PATH = HERE.parent / "knot_census_names.csv"
HEADER_PATH = HERE.parent / "rolfsentable.h"

HYPERBOLIC_DB_PREFIX = "Cusped hyperbolic census"


def bare_census_name(raw: str) -> str:
    """"m004 : #1" -> "m004"; "L104001" -> "L104001" (already bare)."""
    return raw.split(" : ")[0].strip()


def bare_snappy_name(name: str) -> str:
    """"m004(0,0)" -> "m004" (strips SnapPy's Dehn filling suffix)."""
    return name.split("(")[0].strip()


def load_rows():
    rows = []
    with CSV_PATH.open(newline="") as f:
        for name, crossings, db, census_name in csv.reader(f):
            rows.append((name, int(crossings), db, census_name))
    return rows


def main():
    rows = load_rows()

    by_name = {}
    for name, crossings, db, census_name in rows:
        by_name.setdefault(name, []).append((crossings, db, census_name))

    table = {}  # raw census name -> rolfsen name
    collisions = {}  # raw census name -> set of rolfsen names
    mismatches = []  # (rolfsen name, raw census name, snappy identify() names)
    anomalies = []  # rolfsen names with no hit but SnapPy says hyperbolic
    validated = 0

    for name, entries in sorted(by_name.items()):
        hits = [(db, cn) for (_, db, cn) in entries if cn]

        if not hits:
            m = snappy.Manifold(name)
            if m.volume() > 1e-6:
                anomalies.append(name)
            continue

        for db, census_name in hits:
            if db.startswith(HYPERBOLIC_DB_PREFIX):
                m = snappy.Manifold(name)
                ids = {bare_snappy_name(str(x)) for x in m.identify()}
                bare = bare_census_name(census_name)
                if bare not in ids:
                    mismatches.append((name, census_name, sorted(ids)))
                    continue
                validated += 1

            # Keyed by the *base* manifold name (the " : #N" suffix stripped),
            # not the full raw string: entries that share a base name are
            # different ideal triangulations of the same manifold (verified
            # by identical volume across every such group in
            # cusped-hyp-or-census-9), not different manifolds. Matching on
            # the full string would only ever recognize whichever one
            # specific triangulation knotbuilder.cpp + simplify() happened to
            # produce for this PD code, and silently miss the same manifold
            # reached via any other triangulation. Christy hits have no " : "
            # suffix to begin with, so this is a no-op for them.
            #
            # Christy hits have no independent SnapPy cross-check available,
            # so they're trusted as-is; only entries that failed validation
            # above are excluded from the table.
            collisions.setdefault(bare_census_name(census_name), set()).add(name)

    for census_name, names in collisions.items():
        if len(names) > 1:
            continue
        table[census_name] = next(iter(names))

    ambiguous = {k: v for k, v in collisions.items() if len(v) > 1}

    print(f"{len(by_name)} knots, {validated} hyperbolic hits validated against "
          f"SnapPy, {len(mismatches)} mismatches, {len(anomalies)} anomalies, "
          f"{len(ambiguous)} collisions", file=sys.stderr)

    for name, census_name, ids in mismatches:
        print(f"MISMATCH: {name} -> {census_name!r} not confirmed by SnapPy "
              f"(identify() gave {ids})", file=sys.stderr)

    for name in anomalies:
        print(f"ANOMALY: {name} has no census hit, but SnapPy reports it as "
              f"hyperbolic (complement likely needs >9 tetrahedra)",
              file=sys.stderr)

    for census_name, names in ambiguous.items():
        print(f"COLLISION: {census_name!r} maps to multiple Rolfsen names: "
              f"{sorted(names)} -- excluded from the generated table",
              file=sys.stderr)

    write_header(table)
    print(f"wrote {HEADER_PATH} ({len(table)} entries)", file=sys.stderr)

    return 1 if (mismatches or ambiguous) else 0


def write_header(table):
    lines = []
    lines.append("//")
    lines.append("//  rolfsentable.h")
    lines.append("//")
    lines.append("//  GENERATED FILE -- do not edit by hand.")
    lines.append("//  Produced by tests/gen_knot_census_names.cpp + "
                  "tests/gen_knot_names_header.py")
    lines.append("//  from tests/pd_codes_up_to_13_crossings.csv (KnotInfo), "
                  "cross-validated against SnapPy.")
    lines.append("//")
    lines.append("")
    lines.append("#ifndef ROLFSENTABLE_H")
    lines.append("#define ROLFSENTABLE_H")
    lines.append("")
    lines.append("#include <optional>")
    lines.append("#include <string>")
    lines.append("#include <unordered_map>")
    lines.append("")
    lines.append("/*! \\file utils/surfer/rolfsentable.h")
    lines.append(" *  \\brief Translates Regina census hit names (e.g. \"m004 : "
                  "#1\", \"L104001\")")
    lines.append(" *  into classical Rolfsen knot table names (e.g. \"4_1\"), "
                  "for knots up to")
    lines.append(" *  10 crossings -- the largest crossing number Rolfsen's "
                  "table covers.")
    lines.append(" */")
    lines.append("")
    lines.append("namespace rolfsentable {")
    lines.append("")
    lines.append("/**")
    lines.append(" * Base census manifold name (the \" : #N\" suffix stripped) "
                  "-> classical")
    lines.append(" * Rolfsen name. Keyed by base name, not the full raw hit "
                  "string: entries")
    lines.append(" * sharing a base name in Regina's cusped hyperbolic census "
                  "are different")
    lines.append(" * ideal triangulations of the *same* manifold (verified by "
                  "identical")
    lines.append(" * volume across every such group), not different "
                  "manifolds, so matching")
    lines.append(" * on the base name is what makes this table robust to "
                  "which specific")
    lines.append(" * triangulation a caller's simplify() happens to land on.")
    lines.append(" */")
    lines.append("inline const std::unordered_map<std::string, std::string> "
                  "table = {")
    for census_name in sorted(table):
        rolfsen = table[census_name]
        lines.append(f'    {{"{escape(census_name)}", "{escape(rolfsen)}"}},')
    lines.append("};")
    lines.append("")
    lines.append("/**")
    lines.append(" * Returns the Rolfsen name (e.g. \"4_1\") corresponding to "
                  "a raw Regina")
    lines.append(" * census hit name (e.g. \"m004 : #1\" or \"L104001\"), or "
                  "nullopt if unknown.")
    lines.append(" * Any \" : #N\" suffix is stripped before looking up in "
                  "table (see its")
    lines.append(" * comment above), so callers can pass CensusHit::name() "
                  "unmodified.")
    lines.append(" */")
    lines.append("inline std::optional<std::string> "
                  "rolfsenName(const std::string &censusName) {")
    lines.append("    std::string base = censusName.substr(0, "
                  "censusName.find(\" : \"));")
    lines.append("    auto it = table.find(base);")
    lines.append("    if (it == table.end())")
    lines.append("        return std::nullopt;")
    lines.append("    return it->second;")
    lines.append("}")
    lines.append("")
    lines.append("} // namespace rolfsentable")
    lines.append("")
    lines.append("#endif // ROLFSENTABLE_H")
    lines.append("")

    HEADER_PATH.write_text("\n".join(lines))


def escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


if __name__ == "__main__":
    sys.exit(main())

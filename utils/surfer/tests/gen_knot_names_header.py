#!/usr/bin/env python3
"""
Validates ../knot_census_names.csv (produced by gen_knot_census_names.cpp)
and, if present, ../link_census_names.csv (produced by the same tool's
--links mode) against SnapPy, then emits ../linknames.h: a census-name ->
classical-name lookup table (Rolfsen knot names + Thistlethwaite link
names) for identifycomplement.cpp.

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

Knot rows and link rows are validated identically and merged into one
table: SnapPy loads a Thistlethwaite link name (e.g. "L6a1") as a census
manifold exactly the same way it loads a Rolfsen knot name (e.g. "4_1"),
multi-cusp Dehn-filling suffixes and all (confirmed empirically: e.g.
snappy.Manifold("L6a1").identify() includes "s780(0,0)(0,0)", whose bare
name matches this project's own Census::lookup() hit for L6a1 exactly) --
and a non-hyperbolic link (e.g. the Hopf link, "L2a1") loads fine with
near-zero volume, same as a non-hyperbolic Rolfsen name like "3_1" (a torus
knot) already does. So no link-specific validation logic is needed --
knot/link rows only differ in which input CSV and console-report label
they come from. Base census names are trusted to be disjoint between the
two tables in practice (knots and links occupy different census entries);
a genuine collision would still be caught and excluded exactly like a
same-table collision is.

Validation performed:
  - Every row hitting "Cusped hyperbolic census (...)" is cross-checked
    against snappy.Manifold(classical_name).identify(): the bare census
    name (e.g. "m004" from "m004 : #1") must appear among SnapPy's
    independent identification of the same knot/link.
  - Every knot/link with *no* census hit is cross-checked as genuinely
    non-hyperbolic (near-zero volume) via SnapPy. A "no hit" entry that
    SnapPy says *is* hyperbolic is an anomaly -- it just means the
    complement needs more than 9 tetrahedra to fall inside Regina's shipped
    census, not a bug, but it's worth seeing called out explicitly.
  - Collisions (the same base census name mapping to two different
    classical names) are detected and excluded from the generated table
    rather than silently resolved.
"""

import csv
import pathlib
import sys

import snappy

HERE = pathlib.Path(__file__).resolve().parent
KNOT_CSV_PATH = HERE.parent / "knot_census_names.csv"
LINK_CSV_PATH = HERE.parent / "link_census_names.csv"
HEADER_PATH = HERE.parent / "linknames.h"

HYPERBOLIC_DB_PREFIX = "Cusped hyperbolic census"


def bare_census_name(raw: str) -> str:
    """"m004 : #1" -> "m004"; "L104001" -> "L104001" (already bare)."""
    return raw.split(" : ")[0].strip()


def bare_snappy_name(name: str) -> str:
    """"m004(0,0)" -> "m004"; "s780(0,0)(0,0)" -> "s780" (strips SnapPy's
    Dehn filling suffix(es) -- one cusp for a knot, one per cusp for a
    link, but splitting on the first "(" isolates the base name either
    way)."""
    return name.split("(")[0].strip()


def load_rows(csv_path: pathlib.Path):
    """Ignores the trailing isosig column (gen_knot_census_names.cpp's
    complement isoSig, kept for feeding a SnapPy identification pass on
    miss rows -- see tools/identify_link_isosigs.py) -- not needed here."""
    if not csv_path.exists():
        return []
    rows = []
    with csv_path.open(newline="") as f:
        for name, crossings, db, census_name, _isosig in csv.reader(f):
            rows.append((name, int(crossings), db, census_name))
    return rows


def validate(by_name, kind, table, collisions, mismatches, anomalies):
    """Validates one kind's ("knot"/"link") by_name grouping against
    SnapPy, merging confirmed hits into the shared table/collisions dicts.
    Returns the number of hyperbolic hits validated."""
    validated = 0
    for name, entries in sorted(by_name.items()):
        hits = [(db, cn) for (_, db, cn) in entries if cn]

        if not hits:
            m = snappy.Manifold(name)
            if m.volume() > 1e-6:
                anomalies.append((kind, name))
            continue

        for db, census_name in hits:
            if db.startswith(HYPERBOLIC_DB_PREFIX):
                m = snappy.Manifold(name)
                ids = {bare_snappy_name(str(x)) for x in m.identify()}
                bare = bare_census_name(census_name)
                if bare not in ids:
                    mismatches.append((kind, name, census_name, sorted(ids)))
                    continue
                validated += 1

            # Keyed by the *base* manifold name (the " : #N" suffix
            # stripped), not the full raw string -- see this file's own
            # module docstring for why. Christy hits have no " : " suffix
            # to begin with, so this is a no-op for them, and have no
            # independent SnapPy cross-check available, so they're trusted
            # as-is; only entries that failed validation above are
            # excluded from the table.
            collisions.setdefault(bare_census_name(census_name), set()).add(name)
    return validated


def main():
    knot_by_name = {}
    for name, crossings, db, census_name in load_rows(KNOT_CSV_PATH):
        knot_by_name.setdefault(name, []).append((crossings, db, census_name))

    link_by_name = {}
    for name, crossings, db, census_name in load_rows(LINK_CSV_PATH):
        link_by_name.setdefault(name, []).append((crossings, db, census_name))

    table = {}  # raw census name -> classical name
    collisions = {}  # raw census name -> set of classical names
    mismatches = []  # (kind, name, raw census name, snappy identify() names)
    anomalies = []  # (kind, name) with no hit but SnapPy says hyperbolic

    knot_validated = validate(knot_by_name, "knot", table, collisions,
                              mismatches, anomalies)
    link_validated = validate(link_by_name, "link", table, collisions,
                              mismatches, anomalies)

    for census_name, names in collisions.items():
        if len(names) > 1:
            continue
        table[census_name] = next(iter(names))

    ambiguous = {k: v for k, v in collisions.items() if len(v) > 1}

    print(f"{len(knot_by_name)} knots ({knot_validated} hyperbolic hits "
          f"validated), {len(link_by_name)} links ({link_validated} "
          f"hyperbolic hits validated) against SnapPy, {len(mismatches)} "
          f"mismatches, {len(anomalies)} anomalies, {len(ambiguous)} "
          "collisions", file=sys.stderr)

    if not link_by_name:
        print(f"NOTE: {LINK_CSV_PATH} not found -- generating knot-only "
              "(run gen_knot_census_names --links to also cover links)",
              file=sys.stderr)

    for kind, name, census_name, ids in mismatches:
        print(f"MISMATCH: {kind} {name} -> {census_name!r} not confirmed "
              f"by SnapPy (identify() gave {ids})", file=sys.stderr)

    for kind, name in anomalies:
        print(f"ANOMALY: {kind} {name} has no census hit, but SnapPy "
              "reports it as hyperbolic (complement likely needs >9 "
              "tetrahedra)", file=sys.stderr)

    for census_name, names in ambiguous.items():
        print(f"COLLISION: {census_name!r} maps to multiple classical "
              f"names: {sorted(names)} -- excluded from the generated "
              "table", file=sys.stderr)

    write_header(table)
    print(f"wrote {HEADER_PATH} ({len(table)} entries)", file=sys.stderr)

    return 1 if (mismatches or ambiguous) else 0


def write_header(table):
    lines = []
    lines.append("//")
    lines.append("//  linknames.h")
    lines.append("//")
    lines.append("//  GENERATED FILE -- do not edit by hand.")
    lines.append("//  Produced by tests/gen_knot_census_names.cpp + "
                  "tests/gen_knot_names_header.py")
    lines.append("//  from tests/pd_codes_up_to_13_crossings.csv (KnotInfo, "
                  "classical Rolfsen")
    lines.append("//  names) and "
                  "links_4d_smooth_slice_genus_11_crossings_pd_codes.csv")
    lines.append("//  (Thistlethwaite Link Table names), cross-validated "
                  "against SnapPy.")
    lines.append("//")
    lines.append("")
    lines.append("#ifndef LINKNAMES_H")
    lines.append("#define LINKNAMES_H")
    lines.append("")
    lines.append("#include <optional>")
    lines.append("#include <string>")
    lines.append("#include <unordered_map>")
    lines.append("")
    lines.append("/*! \\file utils/surfer/linknames.h")
    lines.append(" *  \\brief Translates Regina census hit names (e.g. \"m004 : "
                  "#1\", \"L104001\")")
    lines.append(" *  into classical names: Rolfsen knot table names (e.g. "
                  "\"4_1\") for knots up")
    lines.append(" *  to 10 crossings -- the largest crossing number "
                  "Rolfsen's table covers --")
    lines.append(" *  and Thistlethwaite Link Table names (e.g. \"L6a1\") "
                  "for links. One shared")
    lines.append(" *  table/namespace rather than a separate one per "
                  "table: knot and link")
    lines.append(" *  census-hit names are disjoint in practice, so "
                  "there's no ambiguity in")
    lines.append(" *  looking both up through the same name().")
    lines.append(" */")
    lines.append("")
    lines.append("namespace linknames {")
    lines.append("")
    lines.append("/**")
    lines.append(" * Base census manifold name (the \" : #N\" suffix stripped) "
                  "-> classical")
    lines.append(" * name (a Rolfsen knot name or a Thistlethwaite link name). "
                  "Keyed by base")
    lines.append(" * name, not the full raw hit string: entries sharing a "
                  "base name in")
    lines.append(" * Regina's cusped hyperbolic census are different ideal "
                  "triangulations of")
    lines.append(" * the *same* manifold (verified by identical volume "
                  "across every such")
    lines.append(" * group), not different manifolds, so matching on the "
                  "base name is what")
    lines.append(" * makes this table robust to which specific "
                  "triangulation a caller's")
    lines.append(" * simplify() happens to land on.")
    lines.append(" */")
    lines.append("inline const std::unordered_map<std::string, std::string> "
                  "table = {")
    for census_name in sorted(table):
        classical = table[census_name]
        lines.append(f'    {{"{escape(census_name)}", "{escape(classical)}"}},')
    lines.append("};")
    lines.append("")
    lines.append("/**")
    lines.append(" * Returns the classical name (a Rolfsen knot name like "
                  "\"4_1\" or a")
    lines.append(" * Thistlethwaite link name like \"L6a1\") corresponding "
                  "to a raw Regina")
    lines.append(" * census hit name (e.g. \"m004 : #1\" or \"L104001\"), or "
                  "nullopt if unknown.")
    lines.append(" * Any \" : #N\" suffix is stripped before looking up in "
                  "table (see its")
    lines.append(" * comment above), so callers can pass CensusHit::name() "
                  "unmodified.")
    lines.append(" */")
    lines.append("inline std::optional<std::string> "
                  "name(const std::string &censusName) {")
    lines.append("    std::string base = censusName.substr(0, "
                  "censusName.find(\" : \"));")
    lines.append("    auto it = table.find(base);")
    lines.append("    if (it == table.end())")
    lines.append("        return std::nullopt;")
    lines.append("    return it->second;")
    lines.append("}")
    lines.append("")
    lines.append("} // namespace linknames")
    lines.append("")
    lines.append("#endif // LINKNAMES_H")
    lines.append("")

    HEADER_PATH.write_text("\n".join(lines))


def escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


if __name__ == "__main__":
    sys.exit(main())

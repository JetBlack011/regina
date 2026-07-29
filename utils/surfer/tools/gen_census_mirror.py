#!/usr/bin/env python3
"""Builds utils/surfer/data/census-mirror.sqlite: a SQLite mirror of the 3
census databases relevant to utils/surfer's boundary-complement lookups
(regina::Census::lookup() on a cusped Triangulation<3> can only ever match
one of these), plus SnapPy-identified isoSigs harvested from
custom_analysis/token_identification_cache.json that Regina's own census
misses.

Queried at runtime by mirrorCensusLookup() (linkcomplement.cpp) via one
read-only SQLite connection per thread -- no hand-rolled mutex needed, unlike
the real regina::Census::lookup(), whose Tokyo Cabinet backend is not safe
under concurrent opens (see linkcomplement.h's censusLookupMutex comment).
A mirror miss falls through to that real, mutex-guarded lookup unchanged, so
this mirror only needs to be *a* correct subset, never the sole source of
truth.

Key format: regina::Census::lookup() queries its on-disk databases by
Triangulation<3>::neoSig() ("second-generation" signature -- see
engine/census/census-impl.h), which is what the .xzsig source files are
actually keyed by despite the generic "isosig" column name. Every other
signature utils/surfer computes or stores -- recognitionCache,
EdgeComplement::identify()'s bare fallback identifier, and therefore
custom_analysis/identify_boundaries.py's token_identification_cache.json --
uses the older, unrelated Triangulation<3>::isoSig() ("first-generation").
So each regina-sourced row is re-keyed here, once, from its native neoSig to
the matching isoSig (via Triangulation3.fromSig(neoSig).isoSig(), using
*this repo's own* just-built Python bindings -- the venv's pip-installed
regina is an older release with no neoSig()/fromSig() support at all).
mirrorCensusLookup() then only ever needs to query by isoSig, matching every
other lookup in this program, with no signature-generation bookkeeping at
runtime. SnapPy-sourced rows need no such conversion: their tokens already
come from EdgeComplement::identify()'s isoSig fallback.

Run by hand, not part of the default build (like tests/gen_knot_census_names
+ tests/gen_knot_names_header.py). Requires the repo to already be built
(needs build/python/regina, for the neoSig->isoSig conversion above -- NOT
utils/surfer/.venv, which lacks it):

    python3 tools/gen_census_mirror.py

Re-run whenever engine/data/census/*.xzsig changes (rare) or after
custom_analysis/identify_boundaries.py extends token_identification_cache.json
with a new search's results -- this script always regenerates the mirror
from scratch, so it stays a pure function of those two checked-in inputs.
"""
import argparse
import json
import lzma
import re
import sqlite3
import sys
from pathlib import Path

SURFER_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SURFER_DIR.parent.parent
DEFAULT_BUILD_PYTHON = REPO_ROOT / "build" / "python"

# The only 3 of Regina's 6 Triangulation<3> census databases that can ever
# match a cusped boundary complement (see census.cpp's CensusCollection<
# Triangulation<3>>::init() -- the other 3 are closed-manifold censuses).
CENSUS_DBS = [
    "cusped-hyp-or-census-9",
    "cusped-hyp-nor-census-9",
    "christy-knots-links",
]

DEFAULT_CENSUS_DIR = REPO_ROOT / "engine" / "data" / "census"
DEFAULT_SNAPPY_CACHE = SURFER_DIR / "custom_analysis" / "token_identification_cache.json"
DEFAULT_OUTPUT = SURFER_DIR / "data" / "census-mirror.sqlite"

# Ported from custom_analysis/identify_boundaries.py: distinguishes an
# already-census-hit token (e.g. "t12656 : #2") from a bare isoSig -- only
# the latter is meaningful as a mirror key.
_HASH_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\s*:\s*#\d+$")


def import_regina(build_python):
    """Imports this repo's own just-built regina Python module (needed for
    the neoSig->isoSig conversion below) -- NOT whatever "regina" a bare
    `import regina` might otherwise resolve to (e.g. utils/surfer/.venv's
    older pip-installed release, which lacks neoSig()/fromSig() entirely)."""
    engine_so = build_python / "regina" / "engine.so"
    if not engine_so.exists():
        sys.exit(
            f"error: {engine_so} not found -- build this repo first "
            f"(this script needs its own freshly-built regina Python "
            f"module, not utils/surfer/.venv's, for neoSig->isoSig "
            f"conversion; see this script's docstring)"
        )
    sys.path.insert(0, str(build_python))
    import regina
    return regina


def load_regina_rows(census_dir, regina):
    """Yields (isosig, name) from each of CENSUS_DBS's .xzsig source file,
    re-keying each entry from its native neoSig (what the file's first
    column, and the real census databases, actually use) to the matching
    isoSig (what everything else in this program uses) -- see this script's
    docstring."""
    for db in CENSUS_DBS:
        path = census_dir / f"{db}.xzsig"
        n = 0
        skipped = 0
        with lzma.open(path, "rt") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                neosig, _, name = line.partition(" ")
                try:
                    isosig = regina.Triangulation3.fromSig(neosig).isoSig()
                except Exception as e:
                    print(f"  warning: {db}: could not convert {neosig!r} "
                          f"({e}), skipping", file=sys.stderr)
                    skipped += 1
                    continue
                yield isosig, name
                n += 1
        print(f"  {db}: {n} entries"
              f"{f' ({skipped} skipped)' if skipped else ''}", file=sys.stderr)


def load_snappy_rows(cache_path):
    """Yields (isosig, name) for cache entries that are genuinely new
    identifications: a bare isoSig token (not "Unknot", not an
    already-census-hit string) that SnapPy (directly or via a bounded
    retriangulation search) resolved to a name."""
    if not cache_path.exists():
        print(f"  no SnapPy cache at {cache_path}, skipping", file=sys.stderr)
        return
    with open(cache_path) as f:
        cache = json.load(f)
    n = 0
    for token, result in cache.items():
        if token == "Unknot" or _HASH_RE.match(token):
            continue
        if result.get("kind") != "identified":
            continue
        name = result.get("best_name")
        if not name:
            continue
        yield token, name
        n += 1
    print(f"  {cache_path.name}: {n} SnapPy-identified isoSigs", file=sys.stderr)


def build(output_path, census_dir, snappy_cache, build_python):
    regina = import_regina(build_python)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    conn = sqlite3.connect(output_path)
    conn.execute(
        "CREATE TABLE census ("
        "  isosig TEXT PRIMARY KEY,"
        "  name TEXT NOT NULL,"
        "  source TEXT NOT NULL"
        ")"
    )

    print("Regina census sources:", file=sys.stderr)
    regina_rows = list(load_regina_rows(census_dir, regina))
    conn.executemany(
        "INSERT OR IGNORE INTO census(isosig, name, source) VALUES (?, ?, 'regina')",
        regina_rows,
    )

    print("SnapPy-extended sources:", file=sys.stderr)
    snappy_rows = list(load_snappy_rows(snappy_cache))
    conn.executemany(
        "INSERT OR IGNORE INTO census(isosig, name, source) VALUES (?, ?, 'snappy')",
        snappy_rows,
    )

    conn.commit()
    total = conn.execute("SELECT COUNT(*) FROM census").fetchone()[0]
    conn.close()

    print(
        f"Wrote {output_path}: {total} rows "
        f"({len(regina_rows)} regina candidates, {len(snappy_rows)} snappy candidates, "
        f"{len(regina_rows) + len(snappy_rows) - total} collisions ignored)",
        file=sys.stderr,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--census-dir", type=Path, default=DEFAULT_CENSUS_DIR)
    ap.add_argument("--snappy-cache", type=Path, default=DEFAULT_SNAPPY_CACHE)
    ap.add_argument("--build-python", type=Path, default=DEFAULT_BUILD_PYTHON,
                     help="path to this repo's built python/ dir "
                          "(containing regina/engine.so), NOT "
                          "utils/surfer/.venv")
    args = ap.parse_args()

    build(args.output, args.census_dir, args.snappy_cache, args.build_python)


if __name__ == "__main__":
    main()

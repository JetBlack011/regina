#!/usr/bin/env python3
"""Runs SnapPy identification (via custom_analysis/identify_boundaries.py's
existing identify_token()/identify_pass() pipeline) over every isoSig in
../link_census_names.csv that Regina's own bundled census missed, merging
results into custom_analysis/token_identification_cache.json -- the same
cache tools/gen_census.py already folds into data/census.sqlite as
source='snappy' rows.

This is deliberately thin: identify_token() already does everything
needed (SnapPy identify(), a Regina-census cross-check, and a bounded
retriangulation search for degenerate SnapPea solutions -- see its own
docstring), and identify_pass() already handles the multiprocessing,
checkpointing, and cache file format. Reusing them here (rather than a
simpler standalone snappy.Manifold(isosig).identify() loop) keeps link
identification exactly as thorough as the existing knot-complement
pipeline, not a weaker one-off.

Usage (run with the venv that has SnapPy installed):
    utils/surfer/.venv/bin/python3 tools/identify_link_isosigs.py

Run after regenerating ../link_census_names.csv (see
tests/gen_knot_census_names.cpp's --links mode) and before
tools/gen_census.py (which reads the cache this writes into and rebuilds
data/census.sqlite from scratch).
"""

import csv
import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
SURFER_DIR = HERE.parent
LINK_CSV_PATH = SURFER_DIR / "link_census_names.csv"
CACHE_PATH = SURFER_DIR / "custom_analysis" / "token_identification_cache.json"

sys.path.insert(0, str(SURFER_DIR / "custom_analysis"))
from identify_boundaries import identify_pass  # noqa: E402


def load_miss_isosigs(csv_path: pathlib.Path) -> set:
    """Reads (name, crossings, db, census_name, isosig) rows -- see
    gen_knot_census_names.cpp's writeRow() -- and returns the isoSigs of
    rows with no Regina census hit (empty db/census_name, still-populated
    isosig)."""
    isosigs = set()
    with csv_path.open(newline="") as f:
        for name, crossings, db, census_name, isosig in csv.reader(f):
            if not census_name and isosig:
                isosigs.add(isosig)
    return isosigs


def main():
    if not LINK_CSV_PATH.exists():
        sys.exit(f"error: {LINK_CSV_PATH} not found -- run "
                 "gen_knot_census_names --links first")

    isosigs = load_miss_isosigs(LINK_CSV_PATH)
    workers = max(1, (os.cpu_count() or 4) - 1)
    print(f"=== identify pass ({len(isosigs)} census-miss link isoSigs, "
          f"{workers} workers) ===", file=sys.stderr)

    cache = identify_pass(isosigs, str(CACHE_PATH), workers)

    kinds = {}
    for token in isosigs:
        r = cache[token]
        kinds[r.kind] = kinds.get(r.kind, 0) + 1
    print(f"link identification summary: {kinds}", file=sys.stderr)


if __name__ == "__main__":
    main()

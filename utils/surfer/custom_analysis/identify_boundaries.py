#!/usr/bin/env python
"""Identify the knot/link complements encoded in the `boundary` column of
13n_65_2_layer_collar.csv, and write a new CSV with the translations added.

Run under the project's .venv (has snappy + regina):
    .venv/bin/python identify_boundaries.py [--limit N] [--workers N]

Pipeline (see /home/john/.claude/plans/you-ll-notice-that-you-re-zany-meadow.md):
  1. scan pass  -- stream the CSV once, collect the set of unique boundary
                   tokens (per-component pieces + whole-link isosigs).
  2. identify pass -- resolve every unique token to a knot/link name (or
                   flag it unresolved), in parallel, with a JSON checkpoint
                   cache so an interrupted run can resume.
  3. write pass -- stream the CSV again, translate each row's boundary
                   field using the cache, write the output CSV.
"""
import argparse
import csv
import json
import multiprocessing
import os
import re
import sys
import time
from dataclasses import dataclass, field, asdict

csv.field_size_limit(sys.maxsize)

INPUT_CSV = "13n_65_2_layer_collar.csv"
OUTPUT_CSV = "13n_65_2_layer_collar_identified.csv"
CACHE_FILE = "token_identification_cache.json"

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

# Top-level boundary-component marker: ", N: " or start-of-string "N: ".
# Note the lack of a space before ':' -- this is what distinguishes it from
# census-hit tokens like "t12656 : #1", which have a space before the colon
# and therefore never match this pattern.
_BOUNDARY_RE = re.compile(r"(?:^|,\s)(\d+):\s(.*?)(?=,\s\d+:\s|$)")
# Trailing whole-link isosig, e.g. "... (pLvLLQLAzQQ...)" at the end of a
# boundary component's content (only present when it has >1 piece).
_WHOLE_RE = re.compile(r"\s\(([^()]+)\)\s*$")
# Census-hit suffix as produced by regina.CensusHit.name(), e.g. "t12656 : #2".
_HASH_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\s*:\s*#\d+$")


def parse_boundary(s):
    """Parse a `boundary` field into {index: {"pieces": [...], "whole": str|None}}."""
    result = {}
    for m in _BOUNDARY_RE.finditer(s):
        idx = int(m.group(1))
        content = m.group(2)
        whole = None
        wm = _WHOLE_RE.search(content)
        if wm:
            whole = wm.group(1)
            content = content[: wm.start()]
        pieces = [p.strip() for p in content.split(", ")]
        result[idx] = {"pieces": pieces, "whole": whole}
    return result


def reassemble_boundary(parsed):
    """Inverse of parse_boundary, for round-trip self-checking."""
    parts = []
    for idx in sorted(parsed):
        entry = parsed[idx]
        content = ", ".join(entry["pieces"])
        if entry["whole"] is not None:
            content += f" ({entry['whole']})"
        parts.append(f"{idx}: {content}")
    return ", ".join(parts)


def strip_census_hash(token):
    m = _HASH_RE.match(token)
    return m.group(1) if m else token


# ---------------------------------------------------------------------------
# Identification
# ---------------------------------------------------------------------------

_ROLFSEN_RE = re.compile(r"^\d+(\^\d+)?_\d+$")
_THISTLETHWAITE_RE = re.compile(r"^[LK]\d+[an]\d+$")


def _pick_best_name(names):
    if not names:
        return None
    for pred in (lambda n: _ROLFSEN_RE.match(n), lambda n: _THISTLETHWAITE_RE.match(n)):
        for n in names:
            if pred(n):
                return n
    return names[0]


@dataclass
class IdentificationResult:
    kind: str  # "unknot" | "identified" | "unresolved"
    names: list = field(default_factory=list)
    best_name: str = None
    method: str = None
    diagnostics: dict = field(default_factory=dict)

    def to_json(self):
        return asdict(self)


def _snappy_identify(manifold):
    """Return (solution_type, names) for a snappy Manifold."""
    st = manifold.solution_type()
    if "degenerate" in st:
        return st, []
    return st, [n.name() for n in manifold.identify()]


def identify_token(token, retriangulate_height=3, retriangulate_budget=8000, retriangulate_seconds=20):
    """Resolve a single boundary token to an IdentificationResult.

    Imports snappy/regina lazily so this function is safe to ship to
    multiprocessing worker processes.
    """
    import snappy
    import regina
    import time as _time

    if token == "Unknot":
        return IdentificationResult(kind="unknot", names=["Unknot"], best_name="Unknot", method="literal")

    base = strip_census_hash(token)

    try:
        m = snappy.Manifold(base)
    except Exception as e:
        return IdentificationResult(
            kind="unresolved", method="snappy_construct_error",
            diagnostics={"error": str(e)[:300]},
        )

    st, names = _snappy_identify(m)
    if names:
        return IdentificationResult(
            kind="identified", names=names, best_name=_pick_best_name(names),
            method="snappy_identify", diagnostics={"solution_type": st},
        )

    if "degenerate" not in st:
        # Positively-oriented (or similar) solution but identify() came up
        # empty -- SnapPy just doesn't recognize it. Try Regina's census as
        # a cross-check before giving up.
        hit_name = _regina_census_name(regina, base)
        if hit_name:
            return IdentificationResult(
                kind="identified", names=[hit_name], best_name=hit_name,
                method="regina_census", diagnostics={"solution_type": st},
            )
        return IdentificationResult(
            kind="unresolved", method="snappy_no_match",
            diagnostics={"solution_type": st, "volume": _safe_volume(m), "homology": str(m.homology())},
        )

    # Degenerate solution -- try a bounded retriangulation search via Regina,
    # re-testing each candidate triangulation with SnapPy. `threads=1` keeps
    # this single-threaded (parallelism already happens one level up, across
    # tokens); retriangulate() has no built-in cap on triangulations explored
    # or time spent, so both a candidate-count budget and a wall-clock budget
    # are enforced in the callback -- bigger triangulations generate
    # candidates much more slowly, so the count alone isn't a reliable bound.
    found = {}
    explored = {"n": 0}
    deadline = _time.time() + retriangulate_seconds

    def _try_candidate(sig, tri):
        explored["n"] += 1
        try:
            cand = snappy.Manifold(sig)
        except Exception:
            cand = None
        if cand is not None:
            cst, cnames = _snappy_identify(cand)
            if cnames:
                found["names"] = cnames
                found["solution_type"] = cst
                return True
        return explored["n"] >= retriangulate_budget or _time.time() >= deadline

    try:
        t = regina.Triangulation3.fromIsoSig(base)
        t.retriangulate(retriangulate_height, 1, _try_candidate)
    except Exception as e:
        return IdentificationResult(
            kind="unresolved", method="retriangulate_error",
            diagnostics={"error": str(e)[:300]},
        )

    if "names" in found:
        return IdentificationResult(
            kind="identified", names=found["names"], best_name=_pick_best_name(found["names"]),
            method="retriangulate", diagnostics={"solution_type": found["solution_type"]},
        )

    hit_name = _regina_census_name(regina, base)
    if hit_name:
        return IdentificationResult(
            kind="identified", names=[hit_name], best_name=hit_name,
            method="regina_census_after_retriangulate", diagnostics={"solution_type": st},
        )

    diag = {"solution_type": st}
    try:
        diag["tetrahedra"] = regina.Triangulation3.fromIsoSig(base).size()
    except Exception:
        pass
    try:
        diag["homology"] = str(m.homology())
    except Exception:
        pass
    return IdentificationResult(kind="unresolved", method="degenerate_unresolved", diagnostics=diag)


def _regina_census_name(regina, base):
    try:
        t = regina.Triangulation3.fromIsoSig(base)
        hits = regina.Census.lookup(t)
        if hits:
            return hits[0].name()
    except Exception:
        pass
    return None


def _safe_volume(m):
    try:
        return float(m.volume())
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def collect_boundary_tokens(row_boundaries):
    """Given an iterable of raw `boundary` strings, return the set of unique
    tokens (pieces + whole-link isosigs) and validate round-tripping."""
    tokens = set()
    for s in row_boundaries:
        parsed = parse_boundary(s)
        rebuilt = reassemble_boundary(parsed)
        if rebuilt != s:
            raise ValueError(f"parser round-trip mismatch:\n  original: {s!r}\n  rebuilt:  {rebuilt!r}")
        for entry in parsed.values():
            for piece in entry["pieces"]:
                tokens.add(piece)
            if entry["whole"] is not None:
                tokens.add(entry["whole"])
    return tokens


def scan_pass(path, limit=None):
    tokens = set()
    n = 0
    t0 = time.time()
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n += 1
            tokens |= collect_boundary_tokens([row["boundary"]])
            if n % 500000 == 0:
                print(f"  scan: {n} rows, {len(tokens)} unique tokens, {time.time()-t0:.1f}s", file=sys.stderr)
            if limit and n >= limit:
                break
    print(f"scan pass done: {n} rows, {len(tokens)} unique tokens, {time.time()-t0:.1f}s", file=sys.stderr)
    return tokens, n


def load_cache(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        raw = json.load(f)
    return {k: IdentificationResult(**v) for k, v in raw.items()}


def save_cache(path, cache):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump({k: v.to_json() for k, v in cache.items()}, f)
    os.replace(tmp, path)


def _worker_identify(token):
    return token, identify_token(token)


def identify_pass(tokens, cache_path, workers, checkpoint_every=200, checkpoint_seconds=15):
    cache = load_cache(cache_path)
    todo = [t for t in tokens if t not in cache]
    print(f"identify pass: {len(cache)} cached, {len(todo)} to compute, {workers} workers", file=sys.stderr)
    if not todo:
        return cache

    t0 = time.time()
    done = 0
    last_checkpoint = t0
    # Small chunksize so a handful of slow (degenerate-retriangulation) tokens
    # in one chunk don't hold back reporting/checkpointing of everything else.
    with multiprocessing.Pool(workers) as pool:
        for token, result in pool.imap_unordered(_worker_identify, todo, chunksize=4):
            cache[token] = result
            done += 1
            now = time.time()
            if done % checkpoint_every == 0 or (now - last_checkpoint) >= checkpoint_seconds:
                save_cache(cache_path, cache)
                last_checkpoint = now
                elapsed = now - t0
                print(f"  identify: {done}/{len(todo)} ({elapsed:.1f}s, {done/elapsed:.1f}/s)", file=sys.stderr)
    save_cache(cache_path, cache)
    print(f"identify pass done: {len(todo)} computed in {time.time()-t0:.1f}s", file=sys.stderr)
    return cache


def _format_piece(token, cache):
    r = cache[token]
    if r.kind == "unknot":
        return "Unknot"
    if r.kind == "identified":
        return r.best_name
    return f"{token} [unresolved]"


def translate_boundary(s, cache):
    parsed = parse_boundary(s)
    all_identified = True
    detail = {}
    out_parts = []
    for idx in sorted(parsed):
        entry = parsed[idx]
        piece_details = []
        piece_names = []
        for piece in entry["pieces"]:
            r = cache[piece]
            if r.kind != "identified" and r.kind != "unknot":
                all_identified = False
            piece_names.append(_format_piece(piece, cache))
            piece_details.append({"token": piece, **r.to_json()})
        content = ", ".join(piece_names)
        whole_detail = None
        if entry["whole"] is not None:
            wr = cache[entry["whole"]]
            if wr.kind != "identified":
                all_identified = False
            whole_detail = {"token": entry["whole"], **wr.to_json()}
            if len(entry["pieces"]) > 1:
                content += f" ({_format_piece(entry['whole'], cache)})"
        out_parts.append(f"{idx}: {content}")
        detail[str(idx)] = {"pieces": piece_details, "whole": whole_detail}
    return ", ".join(out_parts), detail, all_identified


def write_pass(in_path, out_path, cache, limit=None):
    n = 0
    resolved_rows = 0
    t0 = time.time()
    with open(in_path, newline="") as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames + ["boundary_identified", "boundary_identification_json", "all_identified"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            n += 1
            identified_str, detail, all_ident = translate_boundary(row["boundary"], cache)
            row["boundary_identified"] = identified_str
            row["boundary_identification_json"] = json.dumps(detail, separators=(",", ":"))
            row["all_identified"] = all_ident
            if all_ident:
                resolved_rows += 1
            writer.writerow(row)
            if n % 500000 == 0:
                print(f"  write: {n} rows, {time.time()-t0:.1f}s", file=sys.stderr)
            if limit and n >= limit:
                break
    print(f"write pass done: {n} rows written, {resolved_rows} fully identified ({100*resolved_rows/n:.1f}%), {time.time()-t0:.1f}s", file=sys.stderr)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=INPUT_CSV)
    ap.add_argument("--output", default=OUTPUT_CSV)
    ap.add_argument("--cache", default=CACHE_FILE)
    ap.add_argument("--limit", type=int, default=None, help="only process the first N rows (for testing)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    args = ap.parse_args()

    print(f"=== scan pass (limit={args.limit}) ===", file=sys.stderr)
    tokens, n_rows = scan_pass(args.input, limit=args.limit)

    print(f"=== identify pass ({len(tokens)} unique tokens, {args.workers} workers) ===", file=sys.stderr)
    cache = identify_pass(tokens, args.cache, args.workers)

    kinds = {}
    for r in cache.values():
        kinds[r.kind] = kinds.get(r.kind, 0) + 1
    print(f"identification summary: {kinds}", file=sys.stderr)

    print(f"=== write pass -> {args.output} ===", file=sys.stderr)
    write_pass(args.input, args.output, cache, limit=args.limit)


if __name__ == "__main__":
    main()

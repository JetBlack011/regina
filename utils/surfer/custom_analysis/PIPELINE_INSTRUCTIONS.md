# Running the concordance-atlas pipeline on another machine

Four scripts, run in order. Steps 2-3 are the CPU-heavy ones (steps 1 and 4
are plain Python, no snappy/regina needed).

## 0. Files to copy over

From this machine (`/home/john/Projects/triangles/triangle_analysis/`):

- `merge_shards.py`, `identify_boundaries.py`, `build_concordances.py`,
  `export_graph_data.py` — the pipeline scripts
- `graph_viz_template.html` — the visualization template `export_graph_data.py` splices into
- `4d_smooth_slice_genus_13_crossings.csv` — the slice-genus reference table
- Either the raw shard directory (`13n_459_2_layer_collar/`, to merge fresh
  on the new machine — the search is still running here, so the shards keep
  growing) **or** the already-merged snapshot,
  `13n_459_2_layer_collar.csv` (221MB, 1,861,099 rows as of this snapshot)

## 1. Set up the environment

Needs Python 3.14 (or close — that's what this was built/tested on) with
`snappy`, `regina`, and their dependency `cypari` installed. What's in this
project's `.venv`:

```
snappy==3.3.2  regina==7.4.1  cypari==2.5.6  snappy_manifolds==1.4
```

```
python3 -m venv .venv
.venv/bin/pip install snappy regina
```

If `pip install regina` doesn't have a prebuilt wheel for the target
platform, check Regina's own install docs (regina-normal.github.io) for a
platform-specific package/build — that's outside what these scripts can
help with.

## 2. Merge the shards (skip if you copied the already-merged CSV)

The search tool writes one shard file per thread as
`<name>.csv.tmp.<PID>.shard<i>`. A shard directory can contain leftovers
from earlier crashed/restarted runs too — **only merge the PID of the
currently-running process**, or you'll duplicate/corrupt data from stale
runs. Find it with:

```
pgrep -af surfer
```

Then, to sanity-check which PID's shards are actually the live/complete
set (stale ones are usually much smaller):

```
for pid in $(ls 13n_459_2_layer_collar/ | grep -oE '\.tmp\.[0-9]+\.' | tr -d '.' | sed 's/tmp//' | sort -u); do
  echo "PID $pid: $(du -ch 13n_459_2_layer_collar/*.tmp.$pid.shard* 2>/dev/null | tail -1)"
done
```

Then merge:

```
python3 merge_shards.py 13n_459_2_layer_collar --pid <PID> --output 13n_459_2_layer_collar.csv
```

This is defensive about the shards still being actively written to (drops
any shard's last line if it has no trailing newline, since that means it
was caught mid-write) — safe to run against a live search.

## 3. Identify every boundary token (the slow step)

```
.venv/bin/python identify_boundaries.py \
    --input 13n_459_2_layer_collar.csv \
    --output 13n_459_2_layer_collar_identified.csv \
    --workers <N>
```

`--workers` defaults to `nproc - 1`. This dataset has ~1.86M rows (vs. 3M
for the 13n_65 run, which took about 2 hours on 11 workers on an otherwise-
quiet 12-core machine — scale your expectate runtime by core count and
whatever else is competing for CPU on the new machine). Progress is
checkpointed continuously to `token_identification_cache.json`, so it's
safe to interrupt and resume by rerunning the same command.

## 4. Build the concordance index

```
.venv/bin/python build_concordances.py \
    --input 13n_459_2_layer_collar_identified.csv \
    --output concordances_13n459.json
```

Fast (a single pass over the identified CSV, no snappy/regina calls) —
seconds to low minutes.

## 5. Build the visualization

```
python3 export_graph_data.py \
    --input concordances_13n459.json \
    --output graph_data_13n459.json \
    --html-output concordance_atlas_13n459.html \
    --hub K13n459
```

`--hub` is the node id to center the layout on — SnapPy's name for 13n_459
is `K13n459` (same `K{crossings}{a|n}{index}` convention as `K13n65`
before; confirmed live against this data's own isosigs). Plain Python, no
snappy/regina needed for this step. Open `concordance_atlas_13n459.html`
directly in a browser — it's fully self-contained (data inlined, no
external requests), so no server needed.

## Notes

- All four scripts print progress/summary stats to stdout/stderr as they
  run — worth watching for the identification summary (counts of
  identified / unresolved / unknot) and the concordance-count summary.
- If you re-run step 2 later to pick up more shard data (search still
  running), just rerun steps 3-5 after it — `identify_boundaries.py`'s
  cache means only the *new* tokens get computed from scratch.

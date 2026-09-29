# Run scripts (halcyon, 2026-09-28)

The scripts behind the cascade runs of 2026-09-28, kept as they ran. They
assume halcyon's layout: `~/cascade-bench/` holding a built tree `v3/` (or
`v9/`), `census.sqlite`, the pool lists, and `results/`; the atlas at
`~/Projects/cobordism-atlas`.

| script | what |
|---|---|
| `one_target.sh` | one target: `cascadesearch` (constructive, master withheld unless flags add it), then `cascade_check.py` |
| `cascade_v8.sh` | the v8 build and Pool B baseline; defines `pool()`, which the others source |
| `knots_run.sh` | a list in two passes (cascade alone, then with `--master-witnesses`), pausing Pool B's loop meanwhile, then every check |
| `ladder.sh` | a list up three rungs: alone at 50k-surface hops; with the master's witnesses; with them and 1M-surface hops. Then checks and `record_all.sh` |
| `record_all.sh` | every CERTIFIED proof of the named runs into the store (`cascade_record.py`), ready to copy into the atlas's `results/cascade/` |
| `unverified_n.py` | the table knots of N crossings the atlas has not verified, rows or not, as a pool list |
| `make_pool_b.py` | Pool B: known answers the atlas has no constructive bound for, tagged by whether a fixed binary (c4) ever searched them |
| `compare_runs.py`, `phases.py`, `hopdetail.py`, `certnodes.py` | Pool A comparisons, per-hop phase timers, where a hop's search wall goes, a certificate's nodes |
| `build_v9.sh` | a second tree built beside a running one, so a running ladder's binary is never replaced |

Results of these runs: `../../README.md` ("Knots through 8 crossings") and
the atlas's `results/cascade/`.

#!/usr/bin/env python
"""Merge a surfer search's sharded CSV output into one file.

The surfer tool (--threads N) writes N shard files per active run, named
<basename>.csv.tmp.<PID>.shard<i>. A directory can accumulate shards from
multiple past runs (crashes/restarts) as well as the current one -- this
only merges shards matching one specific PID, since mixing runs would
duplicate/corrupt the data.

Defensive about reading a live, still-growing shard: if a shard's last line
has no trailing newline, it's almost certainly mid-write and is dropped
rather than risking a truncated row in the merged output.

Usage:
    python3 merge_shards.py <shard_dir> --pid <PID> --output <merged.csv>

To find the right PID: `pgrep -af surfer` for the live run, or compare
`ls <shard_dir> | grep -oE '\\.tmp\\.[0-9]+\\.' | sort -u` against total size
per PID to spot stale/abandoned ones (`du -ch <shard_dir>/*.tmp.<PID>.shard*`).
"""
import argparse
import glob
import os
import re

HEADER = "orientable,genus,punctures,triangles,condition,boundary\n"


def merge(shard_dir, pid, output_path):
    pattern = os.path.join(shard_dir, f"*.tmp.{pid}.shard*")
    shard_files = glob.glob(pattern)
    if not shard_files:
        raise SystemExit(f"no shards matched {pattern}")

    def shard_num(f):
        return int(re.search(r"shard(\d+)$", f).group(1))
    shard_files.sort(key=shard_num)

    total_lines = 0
    dropped_partial = 0
    with open(output_path, "w") as out:
        out.write(HEADER)
        for f in shard_files:
            with open(f, "rb") as fh:
                data = fh.read()
            text = data.decode("utf-8")
            lines = text.split("\n")
            if lines and lines[-1] == "":
                lines.pop()  # trailing empty string from a clean final newline
            elif not data.endswith(b"\n") and lines:
                dropped_partial += 1
                lines.pop()  # likely mid-write; drop the incomplete tail line
            for line in lines:
                out.write(line + "\n")
            total_lines += len(lines)

    print(f"shards merged: {len(shard_files)}")
    print(f"data rows written: {total_lines}")
    print(f"shards with a dropped (likely mid-write) trailing line: {dropped_partial}")
    print(f"wrote {output_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("shard_dir")
    ap.add_argument("--pid", required=True, help="PID suffix identifying which run's shards to merge")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    merge(args.shard_dir, args.pid, args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Root passes per IDDFS round, from a verifyslicegenus row's progress log.

    round_passes.py <row.err> [...]

The progress block prints "roots exhausted this round: X/R (visits V)",
where V counts root passes cumulatively across rounds. So round k's passes
per root are (V at the end of round k - V at the end of round k-1) / R.

What root_budget_start is calibrated against: a 1M-surface c4 row made about
8 passes per root in round 1 and 10 in round 2 (visits ~1,150 then ~2,500
over 140 roots). With resumed passes each pass is cheaper, but the number of
passes still sets how finely the budget interleaves roots before the surface
target stops the search.
"""
import re
import sys

ESC = re.compile(r'\x1b\[[0-9;]*[A-Za-z]')
ROUND = re.compile(r'iddfs round (\d+)/\d+')
VISITS = re.compile(r'roots exhausted this round: \d+/(\d+) \(visits (\d+)\)')

for path in sys.argv[1:]:
    text = ESC.sub('', open(path, errors='replace').read()).replace('\r', '\n')
    last = {}  # round -> (visits, roots) at its last progress block
    current = None
    for line in text.splitlines():
        m = ROUND.search(line)
        if m:
            current = int(m.group(1))
        m = VISITS.search(line)
        if m and current is not None:
            last[current] = (int(m.group(2)), int(m.group(1)))
    previous = 0
    parts = []
    for k in sorted(last):
        visits, roots = last[k]
        parts.append(f"round {k}: {(visits - previous) / roots:.1f} passes/root")
        previous = visits
    print(f"{path}: " + (", ".join(parts) or "no progress blocks"))

"""Append a deviation row to the DESIGN.md §12 ledger with the real clock time.

Usage: python scripts/ledger_add.py D24 "decision" "why" "bias direction"
Inserts after the last ledger row (a line starting with '| D' or '| C'), stamping the time from the system clock.
"""
import datetime
import sys

ident, decision, why, bias = sys.argv[1:5]
stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
lines = open("record/DESIGN.md").read().splitlines()
last = max(k for k, l in enumerate(lines) if l.startswith(("| D", "| C")) and l.count("|") >= 6)
lines.insert(last + 1, f"| {ident} | {stamp} | {decision} | {why} | {bias} |")
open("record/DESIGN.md", "w").write("\n".join(lines) + "\n")
print(lines[last + 1][:80])

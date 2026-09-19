#!/usr/bin/env python3
"""Numerical regression checks for the bundled, independently reproducible scenarios."""

import csv
import math
import sys
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError("CSV has no data rows")
    return rows


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def verify_inventory(rows: list[dict[str, str]]) -> None:
    require(rows[1]["Pb206_atoms"] == "0", "Bateman cancellation created Pb-206 at the first microsecond")
    for row in rows:
        for name, value in row.items():
            if name.endswith("_atoms"):
                require(math.isfinite(float(value)) and float(value) >= 0.0, f"invalid inventory in {name}")


def verify_activation(rows: list[dict[str, str]]) -> None:
    columns = rows[0].keys()
    require("capture_reactions_per_s" in columns, "capture rate must use reaction-rate units")
    require("capture_rate_bq" not in columns, "legacy Bq capture-rate label remains")
    require("untracked_loss_atoms" in columns, "omitted branch loss must be explicit")
    initial = float(rows[0]["total_atoms"])
    for row in rows:
        reported = float(row["total_atoms"])
        loss = float(row["untracked_loss_atoms"])
        require(abs(reported + loss - initial) / initial <= 5e-12, "activation atom balance drifted")
        if row["phase"] == "cooldown":
            require(float(row["capture_reactions_per_s"]) == 0.0, "capture continued during cooldown")
    require(float(rows[-1]["Pu239_atoms"]) > 1e18, "activation scenario no longer breeds expected Pu-239 inventory")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: verify_outputs.py {inventory|activation} FILE.csv")
    kind, filename = sys.argv[1:]
    rows = read_rows(Path(filename))
    if kind == "inventory":
        verify_inventory(rows)
    elif kind == "activation":
        verify_activation(rows)
    else:
        raise SystemExit(f"unknown test kind: {kind}")


if __name__ == "__main__":
    main()

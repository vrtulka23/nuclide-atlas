#!/usr/bin/env python3
"""Plot a CSV emitted by nuclide-atlas; only numpy and matplotlib are required."""
from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

# The script is meant for CI and headless workstations as well as desktops.
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt


def main(input_name: str, output_name: str) -> None:
    with Path(input_name).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit("inventory CSV is empty")
    # A logarithmic timeline is the point of this plot.  The CSV still carries
    # the exact t=0 row; omit it here because a log axis cannot display zero.
    rows = [row for row in rows if float(row["time_yr"]) > 0.0]
    years = [float(row["time_yr"]) for row in rows]
    atom_columns = [key for key in rows[0] if key.endswith("_atoms")]
    figure, (inventory_ax, activity_ax) = plt.subplots(2, 1, figsize=(11, 8), layout="constrained")
    for column in atom_columns:
        # The closed-form chain has cancellation-level numerical remnants in
        # extremely remote daughters at early times. Do not turn those into
        # misleading horizontal lines on a log plot.
        values = [value if (value := float(row[column])) >= 1.0e6 else math.nan for row in rows]
        inventory_ax.plot(years, values, label=column.removesuffix("_atoms"), linewidth=1.5)
    inventory_ax.set(xlabel="time [Julian years]", ylabel="atoms", xscale="log", yscale="log", title="Decay-chain nuclide inventory")
    inventory_ax.legend(ncol=3, fontsize=8)
    activity_ax.plot(years, [max(float(row["total_activity_bq"]), 1e-300) for row in rows], color="#dc2626", linewidth=2)
    activity_ax.set(xlabel="time [Julian years]", ylabel="activity [Bq]", xscale="log", yscale="log", title="Total activity")
    for axis in (inventory_ax, activity_ax):
        axis.grid(True, which="both", alpha=0.22)
    figure.savefig(output_name, dpi=180)
    print(f"wrote {output_name}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: plot_inventory.py inventory.csv inventory.png")
    main(*sys.argv[1:])

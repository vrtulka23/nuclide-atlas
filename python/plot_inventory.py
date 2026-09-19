#!/usr/bin/env python3
"""Plot a CSV emitted by nuclide-atlas; only numpy and matplotlib are required."""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

# The script is meant for CI and headless workstations as well as desktops.
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np


def _log_time_edges(years: np.ndarray) -> np.ndarray:
    """Return pcolormesh bin edges for strictly positive, log-spaced samples."""
    if len(years) == 1:
        return np.array([years[0] / 1.1, years[0] * 1.1])
    middle = np.sqrt(years[:-1] * years[1:])
    return np.concatenate(([years[0] ** 2 / middle[0]], middle, [years[-1] ** 2 / middle[-1]]))


def _phase_marker(axes: list[plt.Axes], years: np.ndarray, rows: list[dict[str, str]]) -> None:
    """Mark irradiation/cooldown transition without obscuring measurements."""
    if "phase" not in rows[0]:
        return
    transition = next((years[i] for i, row in enumerate(rows) if row["phase"] == "cooldown"), None)
    if transition is None:
        return
    for axis in axes:
        axis.axvline(transition, color="#0f172a", linewidth=1.0, linestyle="--", alpha=0.7)
    axes[0].text(transition, 1.02, "cooldown", transform=axes[0].get_xaxis_transform(), ha="left", va="bottom", fontsize=8)


def main(input_name: str, output_name: str) -> None:
    with Path(input_name).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit("inventory CSV is empty")
    # A logarithmic timeline is the point of this plot.  The CSV still carries
    # the exact t=0 row; omit it here because a log axis cannot display zero.
    rows = [row for row in rows if float(row["time_yr"]) > 0.0]
    years = np.asarray([float(row["time_yr"]) for row in rows])
    atom_columns = [
        key for key in rows[0]
        if key.endswith("_atoms") and key not in {"total_atoms", "untracked_loss_atoms"}
    ]
    has_activity = "total_activity_bq" in rows[0]
    has_q_power = "total_q_power_w" in rows[0]
    capture_column = next(
        (name for name in ("capture_reactions_per_s", "capture_rate_bq") if name in rows[0]), None
    )
    has_capture = capture_column is not None
    has_signals = has_activity or has_q_power
    panel_count = 3 + int(has_signals) + int(has_capture)
    figure, axes_grid = plt.subplots(panel_count, 1, figsize=(13, 2.8 + 2.8 * panel_count), layout="constrained", squeeze=False)
    axes = axes_grid[:, 0]
    labels = [column.removesuffix("_atoms") for column in atom_columns]
    atoms = np.asarray([[max(float(row[column]), 1.0) for row in rows] for column in atom_columns])

    heatmap_ax = axes[0]
    mesh = heatmap_ax.pcolormesh(
        _log_time_edges(years), np.arange(len(labels) + 1), np.log10(atoms), cmap="magma", shading="flat"
    )
    heatmap_ax.set(
        xscale="log",
        xlabel="time [Julian years]",
        ylabel="nuclide",
        title="Inventory evolution — log₁₀(atoms)",
        yticks=np.arange(len(labels)) + 0.5,
        yticklabels=labels,
    )
    heatmap_ax.invert_yaxis()
    colorbar = figure.colorbar(mesh, ax=heatmap_ax, pad=0.01)
    colorbar.set_label("log₁₀ atoms")

    next_axis = 1
    trajectories_ax = axes[next_axis]
    next_axis += 1
    palette = plt.cm.tab20(np.linspace(0.0, 1.0, len(labels)))
    for label, values, color in zip(labels, atoms, palette):
        # Tiny cancellation-level remnants in remote daughters are not physical
        # inventory trajectories and would otherwise dominate the log scale.
        visible = np.where(values >= 1e4, values, np.nan)
        trajectories_ax.plot(years, visible, color=color, linewidth=1.35, label=label)
    trajectories_ax.set(
        xscale="log",
        yscale="log",
        xlabel="time [Julian years]",
        ylabel="atoms",
        title="Nuclide trajectories",
    )
    trajectories_ax.grid(True, which="both", alpha=0.2)
    trajectories_ax.legend(ncol=min(5, max(1, len(labels))), fontsize=7, loc="upper left", frameon=False)

    if has_signals:
        signal_ax = axes[next_axis]
        next_axis += 1
        signal_ax.set(xscale="log", yscale="log", xlabel="time [Julian years]", title="Total decay signals")
        if has_activity:
            signal_ax.plot(years, [max(float(row["total_activity_bq"]), 1e-300) for row in rows], color="#ef4444", linewidth=2.2)
            signal_ax.set_ylabel("activity [Bq]", color="#ef4444")
            signal_ax.tick_params(axis="y", labelcolor="#ef4444")
        if has_q_power:
            power_ax = signal_ax.twinx() if has_activity else signal_ax
            power_ax.plot(years, [max(float(row["total_q_power_w"]), 1e-300) for row in rows], color="#8b5cf6", linewidth=2.0)
            power_ax.set_yscale("log")
            power_ax.set_ylabel("Q-value power [W]", color="#8b5cf6")
            power_ax.tick_params(axis="y", labelcolor="#8b5cf6")
        signal_ax.grid(True, which="both", alpha=0.2)

    if has_capture:
        capture_ax = axes[next_axis]
        next_axis += 1
        capture = np.asarray([float(row[capture_column]) for row in rows])
        active = capture > 0.0
        capture_ax.plot(years[active], capture[active], color="#0891b2", linewidth=2.2)
        capture_ax.fill_between(years[active], capture[active], capture[active].min() * 0.7, color="#67e8f9", alpha=0.32)
        capture_ax.set(xscale="log", yscale="log", xlabel="time [Julian years]", ylabel="captures / s", title="Neutron-capture production rate")
        capture_ax.grid(True, which="both", alpha=0.2)

    composition_ax = axes[next_axis]
    final_atoms = atoms[:, -1]
    order = np.argsort(final_atoms)
    composition_ax.barh(np.asarray(labels)[order], final_atoms[order], color=plt.cm.viridis(np.linspace(0.2, 0.9, len(labels))))
    composition_ax.set(xscale="log", xlabel="atoms at final time", title="Final inventory composition")
    composition_ax.grid(True, which="both", axis="x", alpha=0.2)

    _phase_marker(list(axes[:-1]), years, rows)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    figure.savefig(output_name, dpi=180)
    print(f"wrote {output_name}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: plot_inventory.py inventory.csv inventory.png")
    main(*sys.argv[1:])

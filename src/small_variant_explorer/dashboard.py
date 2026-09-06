"""Four-panel depth, callability, variant-type, and evidence dashboard."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from small_variant_explorer.depth import DepthDataset
from small_variant_explorer.variants import CandidateEvidence, EvidenceThresholds


COLORS = {
    "navy": "#19324d",
    "blue": "#2878b5",
    "teal": "#2a9d8f",
    "amber": "#e9c46a",
    "coral": "#e76f51",
    "gray": "#d9e2ec",
}


def create_dashboard(
    depth_dataset: DepthDataset,
    evidence: Iterable[CandidateEvidence],
    thresholds: EvidenceThresholds,
    output_path: Path,
) -> None:
    candidates = tuple(evidence)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.edgecolor": "#5f6b76",
            "axes.labelcolor": "#263238",
            "xtick.color": "#37474f",
            "ytick.color": "#37474f",
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    figure.patch.set_facecolor("#f7f9fb")
    figure.suptitle(
        "Small-Variant Calling & Evidence QC",
        fontsize=18,
        fontweight="bold",
        color=COLORS["navy"],
    )
    for axis in axes.flat:
        axis.set_facecolor("white")
        axis.grid(axis="y", color="#e8edf2", linewidth=0.8, zorder=0)

    _plot_depth(axes[0, 0], depth_dataset, thresholds)
    _plot_callability(axes[0, 1], depth_dataset)
    _plot_variant_types(axes[1, 0], candidates)
    _plot_evidence(axes[1, 1], candidates, thresholds)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, facecolor=figure.get_facecolor())
    plt.close(figure)


def _plot_depth(axis, dataset: DepthDataset, thresholds: EvidenceThresholds) -> None:
    window_size = 10_000
    centers: list[float] = []
    means: list[float] = []
    offset = 0
    for values in dataset.depths.values():
        for start in range(0, len(values), window_size):
            window = values[start : start + window_size]
            centers.append((offset + start + len(window) / 2) / 1_000_000)
            means.append(sum(window) / len(window))
        offset += len(values)
    axis.plot(centers, means, color=COLORS["blue"], linewidth=1.6)
    axis.axhline(
        thresholds.minimum_depth,
        color=COLORS["amber"],
        linestyle="--",
        linewidth=1.3,
        label=f"Minimum {thresholds.minimum_depth}x",
    )
    axis.axhline(
        thresholds.high_depth_threshold,
        color=COLORS["coral"],
        linestyle=":",
        linewidth=1.3,
        label=f"High-depth boundary {thresholds.high_depth_threshold}x",
    )
    axis.set_title("A. Reliable depth by 10 kb window", loc="left")
    axis.set_xlabel("Reference position (Mb)")
    axis.set_ylabel("Mean depth")
    axis.legend(frameon=False, fontsize=9)


def _plot_callability(axis, dataset: DepthDataset) -> None:
    summary = dataset.summary
    low = summary.positions - summary.positions_ge_minimum
    high = summary.positions_ge_minimum - summary.callable_positions
    values = (summary.callable_positions, low, high)
    labels = ("Callable", "Below minimum", "Above boundary")
    colors = (COLORS["teal"], COLORS["gray"], COLORS["coral"])
    bars = axis.bar(labels, [100 * value / summary.positions for value in values], color=colors)
    axis.bar_label(bars, fmt="%.2f%%", padding=3, fontsize=10)
    axis.set_ylim(0, 105)
    axis.set_ylabel("Reference positions (%)")
    axis.set_title("B. Callable genome", loc="left")
    axis.tick_params(axis="x", rotation=12)


def _plot_variant_types(axis, candidates: tuple[CandidateEvidence, ...]) -> None:
    categories = ("SNP", "insertion", "deletion", "complex")
    accepted = Counter(
        candidate.variant_type for candidate in candidates if candidate.status == "accepted"
    )
    rejected = Counter(
        candidate.variant_type for candidate in candidates if candidate.status == "rejected"
    )
    x = range(len(categories))
    width = 0.38
    axis.bar(
        [value - width / 2 for value in x],
        [accepted[name] for name in categories],
        width,
        label="Accepted",
        color=COLORS["teal"],
    )
    axis.bar(
        [value + width / 2 for value in x],
        [rejected[name] for name in categories],
        width,
        label="Rejected",
        color=COLORS["coral"],
    )
    axis.set_xticks(tuple(x), categories)
    axis.set_ylabel("Normalized candidates")
    axis.set_title("C. Candidate types and disposition", loc="left")
    axis.legend(frameon=False)


def _plot_evidence(
    axis,
    candidates: tuple[CandidateEvidence, ...],
    thresholds: EvidenceThresholds,
) -> None:
    if not candidates:
        axis.text(
            0.5,
            0.5,
            "No normalized candidates",
            ha="center",
            va="center",
            transform=axis.transAxes,
            color="#5f6b76",
            fontsize=13,
        )
    for status, color in (("accepted", COLORS["teal"]), ("rejected", COLORS["coral"])):
        selected = [candidate for candidate in candidates if candidate.status == status]
        if selected:
            axis.scatter(
                [candidate.depth for candidate in selected],
                [100 * candidate.alt_fraction for candidate in selected],
                s=38,
                alpha=0.75,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                label=status.capitalize(),
                zorder=3,
            )
    axis.axhline(
        100 * thresholds.minimum_alt_fraction,
        color=COLORS["amber"],
        linestyle="--",
        linewidth=1.3,
        label=f"Minimum alternate fraction ({100 * thresholds.minimum_alt_fraction:.0f}%)",
    )
    axis.axvline(
        thresholds.minimum_depth,
        color=COLORS["navy"],
        linestyle=":",
        linewidth=1.2,
    )
    axis.set_xlabel("Reliable depth")
    axis.set_ylabel("Alternate fraction (%)")
    axis.set_ylim(-3, 103)
    axis.set_title("D. Per-candidate read evidence", loc="left")
    if candidates:
        axis.legend(frameon=False, fontsize=9)

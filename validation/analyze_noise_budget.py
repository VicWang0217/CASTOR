"""What limits the SNR, and how deep each rig goes with time.

Reads the committed forward atlas `data/lulin_forward_performance.csv` (CASTOR's
forward stack equation over telescope/band/magnitude/frame-count) and turns its
per-term variance columns -- which the summary table collapses to six numbers --
into two things an observer actually reasons about:

  * the noise budget versus source magnitude, so the source -> sky -> read-noise
    crossover is visible, and SLT's correlated background floor is shown as the
    constant chunk it is;
  * the SNR=5 limiting magnitude versus total exposure time for all six
    rig/band combinations, with the 2025-11-06 LOT r' night marked at the
    exposure where the noise model was validated end-to-end.

Run from the repository root::

    uv run --with pandas --with matplotlib python validation/analyze_noise_budget.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "validation/data/lulin_forward_performance.csv"
FIGURE = ROOT / "validation/figures/lulin_noise_budget.png"

BUDGET_FRAMES = 30            # one hour of 120 s frames, for the budget panels
SNR_LIMIT = 5.0
# The end-to-end validation night: 15 x 120 s LOT/SOPHIA r'.
ANCHOR_S = 1800.0
ANCHOR_LABEL = "LOT r' noise model\nvalidated here (obs/pred 1.011)"

TERMS = [
    ("source_variance", "source (shot)", "#f2c14e"),
    ("sky_variance", "sky", "#4a90c0"),
    ("detector_variance", "read + dark", "#7d7d7d"),
    ("flatness_variance", "correlated flatness", "#c0504d"),
]
BAND_COLOR = {"g": "#3a8d3a", "r": "#c0504d", "i": "#8064a2"}
RIG_STYLE = {"LOT": "-", "SLT": "--"}


def budget_panel(ax, df, telescope, band):
    sub = df[(df.telescope == telescope) & (df.band == band)
             & (df.frame_count == BUDGET_FRAMES)].sort_values("ab_magnitude")
    mag = sub["ab_magnitude"].values
    stacks = np.vstack([sub[col].values for col, _, _ in TERMS])
    frac = stacks / stacks.sum(axis=0)
    ax.stackplot(mag, frac, colors=[c for _, _, c in TERMS],
                 labels=[lab for _, lab, _ in TERMS])
    ax.set_xlim(mag.min(), mag.max())
    ax.set_ylim(0, 1)
    ax.set_xlabel("source AB magnitude")
    ax.set_ylabel("fraction of total variance")
    flat = "with 2% flatness floor" if telescope == "SLT" else "no flatness term"
    ax.set_title(f"{telescope} {band}'  noise budget, 1 h stack ({flat})")


def limiting_mag(sub):
    """AB where SNR crosses SNR_LIMIT, per frame count (interp on a mono grid)."""
    out = {}
    for _, g in sub.groupby("frame_count"):
        g = g.sort_values("ab_magnitude")
        snr = g["snr"].values          # descending in magnitude
        mag = g["ab_magnitude"].values  # ascending
        if snr.min() > SNR_LIMIT or snr.max() < SNR_LIMIT:
            continue
        # Interpolate mag(SNR): xp must ascend, so reverse to snr low->high.
        out[g["total_exposure_s"].iloc[0]] = float(
            np.interp(SNR_LIMIT, snr[::-1], mag[::-1]))
    return out


def depth_panel(ax, df):
    anchor_mag = None
    for telescope in ("LOT", "SLT"):
        for band in ("g", "r", "i"):
            sub = df[(df.telescope == telescope) & (df.band == band)]
            curve = limiting_mag(sub)
            if not curve:
                continue
            t = np.array(sorted(curve))
            m = np.array([curve[x] for x in t])
            ax.plot(t, m, RIG_STYLE[telescope], color=BAND_COLOR[band],
                    marker="o", ms=3, lw=1.6, label=f"{telescope} {band}'")
            if telescope == "LOT" and band == "r":
                anchor_mag = float(np.interp(ANCHOR_S, t, m))

    ax.invert_yaxis()  # deeper (fainter) up; do this before annotating
    ax.axvline(ANCHOR_S, color="0.35", ls=":", lw=1.2)
    if anchor_mag is not None:
        ax.plot([ANCHOR_S], [anchor_mag], "*", color="k", ms=13, zorder=5)
        ax.annotate(
            ANCHOR_LABEL, xy=(ANCHOR_S, anchor_mag),
            xytext=(ANCHOR_S * 1.5, anchor_mag - 0.7), fontsize=8.5, va="center",
            arrowprops=dict(arrowstyle="->", color="0.35", lw=1.0),
        )
    ax.set_xscale("log")
    ax.set_xlabel("total exposure time  (s, 120 s frames)")
    ax.set_ylabel(f"SNR={SNR_LIMIT:.0f} limiting AB magnitude")
    ax.set_title("How deep each rig goes  (solid LOT, dashed SLT; deeper = up)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(ncol=2, fontsize=8, framealpha=0.9)


def main() -> None:
    df = pd.read_csv(TABLE)
    fig = plt.figure(figsize=(12.5, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], hspace=0.32, wspace=0.22)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    budget_panel(ax1, df, "LOT", "r")
    budget_panel(ax2, df, "SLT", "r")
    ax2.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8.5,
               title="variance term")
    depth_panel(ax3, df)

    fig.savefig(FIGURE, dpi=140, bbox_inches="tight")
    print(f"wrote {FIGURE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

"""Two per-star views the binned end-to-end result does not show.

Both read the committed publishable table `data/lot_r_endtoend_per_star.csv`
(299 saturation-safe LOT/SOPHIA r' stars, 2025-11-06, fifteen 120 s frames) and
need no observatory FITS.

Panel 1 is the canonical validation scatter: predicted versus observed noise
with the 1:1 line, across three decades of flux.  Panel 2 is the saturation
story stated only in prose elsewhere -- every star's peak pixel against its
aperture flux, with the 16-bit ADC ceiling drawn in.

Run from the repository root::

    uv run --with pandas --with matplotlib python validation/analyze_noise_scatter.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "validation/data/lot_r_endtoend_per_star.csv"
PEAK_TABLE = ROOT / "validation/data/lot_r_endtoend_peak_vs_flux.csv"
FIGURE = ROOT / "validation/figures/lot_r_noise_scatter.png"

# SOPHIA's delivered 16-bit frames cap here, not at the preset full well.
ADC_CEILING_E = 60292.2      # 65535 ADU x 0.92 e-/ADU
PEAK_KEEP_E = 50000.0        # endtoend.py FULL_WELL_KEEP after the saturation fix


def main() -> None:
    df = pd.read_csv(TABLE)
    df["noise_predicted_e"] = df["flux_e"] / df["snr_predicted"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4))

    # --- Panel 1: observed vs predicted noise, 1:1 -----------------------
    etc = df["in_etc_regime"].astype(bool)
    both = np.concatenate([df["noise_predicted_e"].values, df["noise_observed_e"].values])
    lo = np.percentile(both, 0.5) * 0.85
    hi = np.percentile(both, 99.5) * 1.18
    line = np.array([lo, hi])
    ax1.plot(line, line, color="0.35", lw=1.3, zorder=1, label="1:1 (unbiased)")
    ax1.scatter(
        df.loc[~etc, "noise_predicted_e"], df.loc[~etc, "noise_observed_e"],
        s=16, c="#c44", alpha=0.55, edgecolor="none",
        label="above 60 ke- (bright end)", zorder=2,
    )
    ax1.scatter(
        df.loc[etc, "noise_predicted_e"], df.loc[etc, "noise_observed_e"],
        s=16, c="#2a6db0", alpha=0.6, edgecolor="none",
        label="ETC regime (< 60 ke-)", zorder=3,
    )
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)
    ax1.set_aspect("equal")
    ax1.set_xlabel("predicted noise per star  (e-)")
    ax1.set_ylabel("observed noise per star  (e-)")
    ax1.set_title("Noise: observed vs predicted (299 stars)")
    ax1.legend(loc="upper left", fontsize=8, framealpha=0.9)
    med = df.loc[etc, "observed_to_predicted"].median()
    ax1.text(
        0.97, 0.05,
        f"ETC-regime median obs/pred SNR = {med:.3f}",
        transform=ax1.transAxes, ha="right", va="bottom", fontsize=8.5,
        bbox=dict(boxstyle="round", fc="white", ec="0.7"),
    )

    # --- Panel 2: peak pixel vs flux, the saturation ceiling -------------
    # Full star list (incl. saturated), so the ceiling is actually visible.
    pk = pd.read_csv(PEAK_TABLE)
    pinned = pk["n_frames_at_adc_ceiling"] > 0
    dropped = (~pinned) & (pk["peak_max_e"] >= PEAK_KEEP_E)
    kept = ~(pinned | dropped)
    ax2.scatter(
        pk.loc[kept, "flux_e_median"], pk.loc[kept, "peak_max_e"],
        s=16, c="#2a6db0", alpha=0.55, edgecolor="none", label="kept",
    )
    ax2.scatter(
        pk.loc[dropped, "flux_e_median"], pk.loc[dropped, "peak_max_e"],
        s=20, c="#e0900a", alpha=0.8, edgecolor="none",
        label=f"peak >= {PEAK_KEEP_E/1e3:.0f} ke- (cut)",
    )
    ax2.scatter(
        pk.loc[pinned, "flux_e_median"], pk.loc[pinned, "peak_max_e"],
        s=26, c="#c44", alpha=0.85, edgecolor="none",
        label="pinned at ADC ceiling",
    )
    ax2.axhline(ADC_CEILING_E, color="k", ls="--", lw=1.2)
    ax2.axhline(PEAK_KEEP_E, color="#e0900a", ls=":", lw=1.1)
    ax2.text(
        pk["flux_e_median"].min() * 1.2, ADC_CEILING_E * 1.04,
        "16-bit ADC ceiling  60,292 e-", fontsize=8.5, va="bottom",
    )
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("aperture flux  (e-)")
    ax2.set_ylabel("brightest pixel across 15 frames  (e-)")
    ax2.set_title("Why the bright end breaks: the ADC ceiling")
    ax2.legend(loc="lower right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(FIGURE, dpi=140)
    print(f"wrote {FIGURE.relative_to(ROOT)}")
    print(f"{int(pinned.sum())} stars pinned at the ADC ceiling")


if __name__ == "__main__":
    main()

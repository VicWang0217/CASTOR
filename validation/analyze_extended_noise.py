"""Extended-source noise versus aperture radius, model against measurement.

The NGC 3621 write-up quotes SNR ratios in a six-row table; this draws the noise
itself.  The measured aperture noise comes from the model-free half-stack
difference (SLT r', 42 x 30 s, 2024-04-12).  Two model curves are shown: the
noise CASTOR predicted when the test was run (Andor read noise 3.3 e-, no
correlated term) and the noise it predicts after the fix on branch
`fix/slt-noise-model` (read noise 9.28 e- plus the 2% background-flatness
floor).  The shipped curve grows only as sqrt(area) and falls away; the fixed
curve tracks the measurement to ~1% once the flatness term dominates past ~4".

Values are the committed per-radius aggregate
`data/extended_ngc3621_noise_vs_radius.csv`; no observatory FITS is needed.

Run from the repository root::

    uv run --with pandas --with matplotlib python validation/analyze_extended_noise.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "validation/data/extended_ngc3621_noise_vs_radius.csv"
FIGURE = ROOT / "validation/figures/extended_noise_vs_radius.png"


def main() -> None:
    df = pd.read_csv(TABLE).sort_values("r_arcsec")
    r = df["r_arcsec"].values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))

    # --- Panel 1: noise vs radius ---------------------------------------
    ax1.plot(r, df["sigma_meas_e"], "o-", color="k", lw=1.8, ms=6,
             label="measured (half-stack difference)")
    ax1.plot(r, df["sigma_pred_shipped_e"], "s--", color="#c44", lw=1.6, ms=5,
             label="shipped model (RN 3.3, no flatness)")
    ax1.plot(r, df["sigma_pred_fixed_e"], "^-", color="#2a6db0", lw=1.6, ms=5,
             label="fixed model (RN 9.28 + 2% flatness)")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("aperture radius  (arcsec)")
    ax1.set_ylabel("noise per aperture, one 30 s frame  (e-)")
    ax1.set_title("Extended-source noise vs aperture (NGC 3621, SLT r')")
    ax1.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
    ax1.grid(True, which="both", alpha=0.25)

    # --- Panel 2: measured/predicted ratio ------------------------------
    ax2.axhline(1.0, color="0.4", lw=1.2)
    ax2.plot(r, df["sigma_meas_e"] / df["sigma_pred_shipped_e"], "s--",
             color="#c44", lw=1.6, ms=5, label="measured / shipped")
    ax2.plot(r, df["sigma_meas_e"] / df["sigma_pred_fixed_e"], "^-",
             color="#2a6db0", lw=1.6, ms=5, label="measured / fixed")
    ax2.set_xscale("log")
    ax2.set_xlabel("aperture radius  (arcsec)")
    ax2.set_ylabel("measured noise / predicted noise")
    ax2.set_title("The shipped model runs away; the fix tracks to ~1%")
    ax2.set_ylim(0.8, 3.4)
    ax2.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
    ax2.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    fig.savefig(FIGURE, dpi=140)
    print(f"wrote {FIGURE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

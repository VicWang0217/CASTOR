"""Measure the LOT/SOPHIA r' colour term in the end-to-end field.

This uses the publishable per-star table made by
``analyze_endtoend_residuals.py``.  A Pan-STARRS PSF-minus-Kron cut removes
extended sources before fitting the instrumental zero point against g-r.

Run from the repository root::

    uv run --with pandas --with matplotlib python validation/analyze_lot_color_term.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lulin


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "validation/data/lot_r_endtoend_per_star.csv"
FIGURE = ROOT / "validation/figures/lot_r_color_term.png"
REPORT = ROOT / "validation/report_sections/color_term.md"

MIN_SNR = 20.0
PSF_KRON_RANGE = (-0.20, 0.05)


def robust_line(x: np.ndarray, y: np.ndarray, reference: float):
    """Fit y = intercept + slope * (x-reference), clipping at three MAD."""
    keep = np.ones(len(x), dtype=bool)
    for _ in range(10):
        design = np.column_stack([np.ones(keep.sum()), x[keep] - reference])
        coefficients = np.linalg.lstsq(design, y[keep], rcond=None)[0]
        residual = y - (coefficients[0] + coefficients[1] * (x - reference))
        centre = np.median(residual[keep])
        scatter = 1.4826 * np.median(np.abs(residual[keep] - centre))
        new_keep = np.abs(residual - centre) < 3.0 * scatter
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return coefficients, keep, scatter


def bootstrap_fit(x: np.ndarray, y: np.ndarray, reference: float, draws: int = 5000):
    rng = np.random.default_rng(20260911)
    estimates = np.empty((draws, 2))
    for index in range(draws):
        sample = rng.integers(0, len(x), len(x))
        estimates[index] = robust_line(x[sample], y[sample], reference)[0]
    return estimates


def analyse():
    table = pd.read_csv(TABLE)
    table["instrumental_zp"] = (
        table["ps1_r_mag"]
        + 2.5 * np.log10(table["flux_e"] / table["exposure_seconds"])
    )
    table["psf_minus_kron"] = table["ps1_r_mag"] - table["ps1_r_kron_mag"]
    candidates = table[
        table["ps1_g_minus_r"].notna()
        & (table["snr_observed"] >= MIN_SNR)
        & (table["ps1_r_kron_mag"] > 0)
        & (table["ps1_r_kron_mag_error"] > 0)
        & table["psf_minus_kron"].between(*PSF_KRON_RANGE)
    ].copy()

    colour = candidates["ps1_g_minus_r"].to_numpy()
    zero_point = candidates["instrumental_zp"].to_numpy()
    reference = float(np.median(colour))
    coefficients, keep, _ = robust_line(colour, zero_point, reference)
    accepted = candidates.iloc[np.flatnonzero(keep)].copy()
    bootstrap = bootstrap_fit(colour, zero_point, reference)
    ci = np.percentile(bootstrap, [2.5, 97.5], axis=0)

    accepted["colour_corrected_residual"] = (
        accepted["instrumental_zp"]
        - coefficients[0]
        - coefficients[1] * (accepted["ps1_g_minus_r"] - reference)
    )
    residual_scatter = float(accepted["colour_corrected_residual"].std(ddof=2))
    q05, q95 = np.quantile(accepted["ps1_g_minus_r"], [0.05, 0.95])
    colour_span_mag = float(coefficients[1] * (q95 - q05))
    colour_span_flux = float(10 ** (0.4 * colour_span_mag) - 1)
    airmass = float(accepted["airmass"].median())
    multi_night_expected = lulin.MEASURED["r"]["zp0"] - lulin.MEASURED["r"]["k"] * airmass

    result = {
        "candidates": len(candidates),
        "accepted": len(accepted),
        "reference": reference,
        "intercept": float(coefficients[0]),
        "intercept_ci": tuple(ci[:, 0]),
        "slope": float(coefficients[1]),
        "slope_ci": tuple(ci[:, 1]),
        "scatter": residual_scatter,
        "q05": float(q05),
        "q95": float(q95),
        "colour_span_mag": colour_span_mag,
        "colour_span_flux": colour_span_flux,
        "multi_night_expected": float(multi_night_expected),
        "airmass": airmass,
    }
    return table, candidates, accepted, bootstrap, result


def make_figure(table, candidates, accepted, bootstrap, result):
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    accent = "#2563eb"
    muted = "#94a3b8"

    ax = axes[0]
    colour_quality = table[table["ps1_g_minus_r"].notna() & (table["snr_observed"] >= MIN_SNR)]
    ax.scatter(colour_quality["psf_minus_kron"], colour_quality["instrumental_zp"],
               s=20, alpha=0.45, color=muted, edgecolors="none")
    ax.axvspan(*PSF_KRON_RANGE, color=accent, alpha=0.10, label="point-source selection")
    ax.axvline(PSF_KRON_RANGE[1], color=accent, linewidth=1)
    ax.set(xlabel="Pan-STARRS r PSF − Kron (mag)", ylabel="Per-star zero point (mag)",
           title="A  Extended sources form the high-ZP branch", xlim=(-0.25, 1.35))
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    ax.scatter(candidates["ps1_g_minus_r"], candidates["instrumental_zp"],
               s=22, alpha=0.30, color=muted, edgecolors="none", label="candidates")
    ax.scatter(accepted["ps1_g_minus_r"], accepted["instrumental_zp"],
               s=22, alpha=0.62, color=accent, edgecolors="none", label="3-MAD fit sample")
    grid = np.linspace(accepted["ps1_g_minus_r"].min(), accepted["ps1_g_minus_r"].max(), 200)
    lines = bootstrap[:, 0, None] + bootstrap[:, 1, None] * (grid[None, :] - result["reference"])
    ax.fill_between(grid, *np.percentile(lines, [2.5, 97.5], axis=0), color=accent, alpha=0.16)
    ax.plot(grid, result["intercept"] + result["slope"] * (grid - result["reference"]),
            color=accent, linewidth=2)
    ax.set(xlabel="Pan-STARRS g−r (mag)", ylabel="Per-star zero point (mag)",
           title="B  LOT r' natural-system colour term")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[2]
    ax.scatter(accepted["ps1_r_mag"], accepted["colour_corrected_residual"],
               s=22, alpha=0.48, color=accent, edgecolors="none")
    accepted = accepted.copy()
    accepted["mag_bin"] = pd.qcut(accepted["ps1_r_mag"], 6, duplicates="drop")
    grouped = accepted.groupby("mag_bin", observed=True)
    ax.plot(grouped["ps1_r_mag"].median(), grouped["colour_corrected_residual"].median(),
            "o-", color="#0f172a", linewidth=1.5, label="magnitude-bin medians")
    ax.axhline(0, color="#0f172a", linewidth=1, linestyle="--")
    ax.set(xlabel="Pan-STARRS r PSF magnitude", ylabel="Colour-corrected ZP residual (mag)",
           title="C  No remaining magnitude-scale trend")
    ax.legend(frameon=False, loc="lower left")

    fig.suptitle("LOT / SOPHIA / r' photometric colour-term audit", fontsize=15, fontweight="bold")
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=180, facecolor="white")
    plt.close(fig)


def write_report(result):
    slope_percent = 100 * (10 ** (0.4 * result["slope"]) - 1)
    zp_difference = result["intercept"] - result["multi_night_expected"]
    text = f"""# LOT r' photometric colour-term audit

The end-to-end field also tests whether calibrating LOT's natural Astrodon r'
response directly against Pan-STARRS r introduces a stellar-colour bias. This
is a count-rate calibration test, separate from the SNR residual test in the
per-star end-to-end section.

![LOT r colour-term analysis](figures/lot_r_color_term.png)

## Selection

The fit starts from the 299 saturation-safe stars, requires measured SNR >=
{MIN_SNR:.0f}, valid Pan-STARRS g/r/i colours, and uses
`-0.20 < r_PSF - r_Kron < 0.05` to select point sources. That last cut matters:
extended objects measured with a circular aperture form a separate zero-point
branch up to 1.4 mag high. Of {result['candidates']} point-source candidates,
{result['accepted']} survive a three-MAD residual clip.

## Result

At the sample's median colour, `g-r = {result['reference']:.3f}`, the observed
zero point is:

`ZP_r = {result['intercept']:.4f} + ({result['slope']:+.4f}) * ((g-r) - {result['reference']:.3f})`

| Quantity | Estimate | Bootstrap 95% CI |
|---|---:|---:|
| Zero point at reference colour | {result['intercept']:.4f} mag | {result['intercept_ci'][0]:.4f} to {result['intercept_ci'][1]:.4f} |
| Colour coefficient | **{result['slope']:+.4f} mag per mag** | {result['slope_ci'][0]:+.4f} to {result['slope_ci'][1]:+.4f} |
| Residual scatter | {result['scatter']:.4f} mag | — |

The coefficient is small but resolved in this sample. One magnitude of g-r
changes the inferred count rate by {slope_percent:+.2f}%. Across the accepted
sample's central 90% colour range ({result['q05']:.2f} to {result['q95']:.2f}),
the full effect is {result['colour_span_mag']:.3f} mag, or
{100 * result['colour_span_flux']:.2f}% in count rate.

The multi-night calibration in `lulin.py` predicts ZP =
{result['multi_night_expected']:.4f} at this sample's median airmass
({result['airmass']:.3f}); this fit differs by {zp_difference:+.4f} mag. That is
an internal consistency check, not independent evidence, because both use the
same LOT field.

## Consequence

This is too small to explain LOT's factor-1.8 r' throughput excess over g' and
i'. It does show that future throughput fits should include a colour term and a
point-source morphology cut. The present r' throughput need not move: its
reference-colour zero point agrees with the existing multi-night relation to
well below one percent. A multi-band fit over independent fields is still
needed before adopting coefficients in the calculator itself.

Regenerate with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_lot_color_term.py
```
"""
    REPORT.write_text(text, encoding="utf-8")


def main():
    table, candidates, accepted, bootstrap, result = analyse()
    make_figure(table, candidates, accepted, bootstrap, result)
    write_report(result)
    print(f"colour term {result['slope']:+.4f} mag/mag ")
    print(f"95% CI {result['slope_ci'][0]:+.4f} to {result['slope_ci'][1]:+.4f}")
    print(f"wrote {FIGURE.relative_to(ROOT)} and {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

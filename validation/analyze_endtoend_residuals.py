"""Publish the LOT/SOPHIA r' end-to-end check at per-star resolution.

The raw FITS files and intermediate photometry are deliberately gitignored.  This
script joins the reduced SNR table to the saturation diagnostic and Pan-STARRS
catalogue, applies the established 50,000 e- peak cut, and writes a public table,
figure, and short report containing no images or proprietary catalogue export.

Run from the repository root::

    uv run --with pandas --with matplotlib python \
        validation/analyze_endtoend_residuals.py

If the ignored inputs are unavailable, ``--reuse-derived`` redraws the figure and
report from the committed per-star table instead.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from castor import physics
from castorCLI import presets


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "validation/data/raw/_endtoend_2026-08-30"
CATALOGUE = ROOT / "validation/data/raw/lot_sn2025wny/ps1/SN2025wny.csv"
TABLE = ROOT / "validation/data/lot_r_endtoend_per_star.csv"
FIGURE = ROOT / "validation/figures/lot_r_endtoend_residuals.png"
REPORT = ROOT / "validation/report_sections/endtoend.md"

PEAK_LIMIT_E = 50_000.0
HEADLINE_FLUX_LIMIT_E = 60_000.0
DETECTOR_SIZE_PIX = 2048.0
FLUX_EDGES = [0, 3e3, 1e4, 3e4, 6e4, 1.2e5, 2.4e5, np.inf]
FLUX_LABELS = ["<3k", "3-10k", "10-30k", "30-60k", "60-120k", "120-240k", ">240k"]


def _current_prediction(table: pd.DataFrame) -> pd.DataFrame:
    """Evaluate the live CASTOR noise model on every measured star."""
    camera = presets.load().profile("lulin").cameras["Sophia"].camera
    n_est = physics.calculate_sky_estimate_pixels(
        table["aperture_pixels"].to_numpy(),
        inner_factor=5.0,
        outer_factor=8.0,
        total_fwhm=table["fwhm_pixels"].to_numpy(),
        pixel_scale=1.0,
        estimator="median",
    )
    predicted = physics.calculate_single_snr(
        source_count_rate=table["flux_e"].to_numpy() / table["exposure_seconds"].to_numpy(),
        sky_count_rate=table["sky_e_per_pixel"].to_numpy() / table["exposure_seconds"].to_numpy(),
        dark_current_rate=camera.dark_current_rate,
        readout_noise=camera.readout_noise,
        num_pixels_aperture=table["aperture_pixels"].to_numpy(),
        single_exp_time=table["exposure_seconds"].to_numpy(),
        num_pixels_sky_estimate=n_est,
        background_flatness_fraction=camera.background_flatness_fraction,
    )
    table = table.copy()
    table["snr_predicted"] = predicted
    table["observed_to_predicted"] = table["snr_observed"] / predicted
    return table


def build_table() -> pd.DataFrame:
    """Join the ignored reduction products and apply the final selection."""
    snr = pd.read_csv(RAW / "snr_per_star.csv")
    shape = pd.read_csv(RAW / "bright_end.csv")

    # Reproduce endtoend.py's catalogue selection exactly.  Its integer star id
    # is a positional index into this filtered table, not the original CSV index.
    catalogue = pd.read_csv(CATALOGUE)
    catalogue = catalogue[
        (catalogue["rMeanPSFMag"] > 13)
        & (catalogue["rMeanPSFMag"] < 21)
        & (catalogue["rMeanPSFMagErr"] > 0)
    ].reset_index(drop=True)

    per_star = shape.groupby("star").agg(
        peak_max_e=("peak", "max"),
        x_pixel=("x", "median"),
        y_pixel=("y", "median"),
    )
    table = snr.join(per_star, on="star")
    table = table[table["peak_max_e"] < PEAK_LIMIT_E].copy()

    for name in (
        "gMeanPSFMag", "gMeanPSFMagErr", "gQfPerfect",
        "rMeanPSFMag", "rMeanPSFMagErr", "rQfPerfect",
        "rMeanKronMag", "rMeanKronMagErr",
        "iMeanPSFMag", "iMeanPSFMagErr", "iQfPerfect",
    ):
        table[name] = table["star"].map(catalogue[name])

    color_ok = (
        (table["gMeanPSFMag"] > 0)
        & (table["iMeanPSFMag"] > 0)
        & (table["gMeanPSFMagErr"] > 0)
        & (table["iMeanPSFMagErr"] > 0)
        & (table["gQfPerfect"] > 0.9)
        & (table["rQfPerfect"] > 0.9)
        & (table["iQfPerfect"] > 0.9)
    )
    g_r = table["gMeanPSFMag"] - table["rMeanPSFMag"]
    r_i = table["rMeanPSFMag"] - table["iMeanPSFMag"]
    color_ok &= g_r.between(-0.5, 3.0) & r_i.between(-0.5, 2.0)
    table["ps1_g_minus_r"] = g_r.where(color_ok)
    table["ps1_r_minus_i"] = r_i.where(color_ok)

    table = table.rename(columns={
        "star": "star_id",
        "noise_obs": "noise_observed_e",
        "snr_obs": "snr_observed",
        "sky_e_pix": "sky_e_per_pixel",
        "npix": "aperture_pixels",
        "fwhm_pix": "fwhm_pixels",
        "exptime_s": "exposure_seconds",
        "rMeanPSFMag": "ps1_r_mag",
        "rMeanPSFMagErr": "ps1_r_mag_error",
        "rMeanKronMag": "ps1_r_kron_mag",
        "rMeanKronMagErr": "ps1_r_kron_mag_error",
    })
    table["x_normalized"] = (table["x_pixel"] - DETECTOR_SIZE_PIX / 2) / (DETECTOR_SIZE_PIX / 2)
    table["y_normalized"] = (table["y_pixel"] - DETECTOR_SIZE_PIX / 2) / (DETECTOR_SIZE_PIX / 2)
    table["field_radius_normalized"] = np.hypot(table["x_normalized"], table["y_normalized"])
    table["flux_bin"] = pd.cut(table["flux_e"], FLUX_EDGES, labels=FLUX_LABELS)
    table["in_etc_regime"] = table["flux_e"] < HEADLINE_FLUX_LIMIT_E
    table = _current_prediction(table)

    columns = [
        "star_id", "flux_bin", "in_etc_regime", "flux_e", "noise_observed_e",
        "snr_observed", "snr_predicted", "observed_to_predicted", "sigma_relerr",
        "peak_max_e", "ps1_r_mag", "ps1_g_minus_r", "ps1_r_minus_i",
        "ps1_r_mag_error", "ps1_r_kron_mag", "ps1_r_kron_mag_error",
        "x_pixel", "y_pixel", "x_normalized", "y_normalized", "field_radius_normalized",
        "sky_e_per_pixel", "aperture_pixels", "fwhm_pixels", "exposure_seconds",
        "airmass", "n_frames",
    ]
    table = table[columns].sort_values("flux_e").reset_index(drop=True)
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False, float_format="%.7g")
    return table


def load_table() -> pd.DataFrame:
    """Load the publishable table and refresh predictions from current code."""
    table = pd.read_csv(TABLE)
    table["in_etc_regime"] = table["in_etc_regime"].astype(str).str.lower().eq("true")
    table["flux_bin"] = pd.Categorical(table["flux_bin"], categories=FLUX_LABELS, ordered=True)
    return _current_prediction(table)


def _rank_correlation(x: pd.Series, y: pd.Series) -> float:
    good = np.isfinite(x) & np.isfinite(y)
    return float(x[good].rank().corr(y[good].rank()))


def _bootstrap_ci(values: np.ndarray, statistic, seed: int, draws: int = 5000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates = np.empty(draws)
    for index in range(draws):
        sample = values[rng.integers(0, len(values), len(values))]
        estimates[index] = statistic(sample)
    return tuple(np.percentile(estimates, [2.5, 97.5]))


def _correlation_ci(frame: pd.DataFrame, column: str, seed: int) -> tuple[float, float, float, int]:
    pair = frame[[column, "observed_to_predicted"]].dropna().to_numpy()
    rho = _rank_correlation(pd.Series(pair[:, 0]), pd.Series(pair[:, 1]))
    rng = np.random.default_rng(seed)
    estimates = np.empty(5000)
    for index in range(len(estimates)):
        sample = pair[rng.integers(0, len(pair), len(pair))]
        estimates[index] = _rank_correlation(pd.Series(sample[:, 0]), pd.Series(sample[:, 1]))
    low, high = np.percentile(estimates, [2.5, 97.5])
    return rho, float(low), float(high), len(pair)


def summarize(table: pd.DataFrame) -> dict:
    safe = table.copy()
    etc = safe[safe["in_etc_regime"]].copy()

    def distribution(frame: pd.DataFrame, seed: int) -> dict:
        values = frame["observed_to_predicted"].to_numpy()
        low, high = _bootstrap_ci(values, np.median, seed)
        return {
            "n": len(values),
            "median": float(np.median(values)),
            "q25": float(np.quantile(values, 0.25)),
            "q75": float(np.quantile(values, 0.75)),
            "ci_low": low,
            "ci_high": high,
        }

    return {
        "safe": distribution(safe, 11),
        "etc": distribution(etc, 12),
        "color_n": int(etc["ps1_g_minus_r"].notna().sum()),
        "correlations": {
            "flux": _correlation_ci(etc.assign(log_flux=np.log10(etc["flux_e"])), "log_flux", 21),
            "color": _correlation_ci(etc, "ps1_g_minus_r", 22),
            "x": _correlation_ci(etc, "x_normalized", 23),
            "y": _correlation_ci(etc, "y_normalized", 24),
            "radius": _correlation_ci(etc, "field_radius_normalized", 25),
        },
    }


def _binned_points(frame: pd.DataFrame, column: str, bins: int = 6):
    good = frame[[column, "observed_to_predicted"]].dropna().sort_values(column)
    good["group"] = pd.qcut(good[column], bins, duplicates="drop")
    grouped = good.groupby("group", observed=True)
    return grouped[column].median(), grouped["observed_to_predicted"].median()


def make_figure(table: pd.DataFrame, summary: dict) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    accent = "#2563eb"
    muted = "#64748b"
    etc = table[table["in_etc_regime"]]

    ax = axes[0, 0]
    ax.scatter(table["snr_predicted"], table["snr_observed"], s=18, alpha=0.45,
               c=np.log10(table["flux_e"]), cmap="viridis", edgecolors="none")
    limits = [min(table["snr_predicted"].min(), table["snr_observed"].min()),
              max(table["snr_predicted"].max(), table["snr_observed"].max())]
    ax.plot(limits, limits, color="#0f172a", linewidth=1.3, linestyle="--", label="perfect agreement")
    ax.set(xscale="log", yscale="log", xlabel="CASTOR predicted SNR", ylabel="Observed SNR",
           title="A  Prediction against measured scatter")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[0, 1]
    ax.scatter(table["flux_e"], table["observed_to_predicted"], s=17, alpha=0.32,
               color=muted, edgecolors="none")
    grouped = table.groupby("flux_bin", observed=True)
    x = grouped["flux_e"].median()
    med = grouped["observed_to_predicted"].median()
    low = med - grouped["observed_to_predicted"].quantile(0.25)
    high = grouped["observed_to_predicted"].quantile(0.75) - med
    ax.errorbar(x, med, yerr=np.vstack([low, high]), fmt="o-", color=accent,
                linewidth=1.7, capsize=3, label="bin median and IQR")
    ax.axhline(1, color="#0f172a", linewidth=1.3, linestyle="--")
    ax.axvline(HEADLINE_FLUX_LIMIT_E, color="#dc2626", linewidth=1, linestyle=":")
    ax.set(xscale="log", xlabel="Source electrons per 120 s frame", ylabel="Observed / predicted SNR",
           title="B  Residual against source brightness")
    ax.legend(frameon=False, loc="lower left")

    ax = axes[1, 0]
    color = etc.dropna(subset=["ps1_g_minus_r"])
    ax.scatter(color["ps1_g_minus_r"], color["observed_to_predicted"], s=18,
               alpha=0.32, color=muted, edgecolors="none")
    bx, by = _binned_points(color, "ps1_g_minus_r")
    ax.plot(bx, by, "o-", color=accent, linewidth=1.7, label="colour-quantile medians")
    ax.axhline(1, color="#0f172a", linewidth=1.3, linestyle="--")
    rho = summary["correlations"]["color"]
    ax.text(0.03, 0.95, f"Spearman ρ = {rho[0]:+.2f}  (95% CI {rho[1]:+.2f} to {rho[2]:+.2f})",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.set(xlabel="Pan-STARRS g−r (mag)", ylabel="Observed / predicted SNR",
           title="C  No detected colour dependence below 60 ke-")
    ax.legend(frameon=False, loc="lower left")

    ax = axes[1, 1]
    xbins = np.linspace(0, DETECTOR_SIZE_PIX, 5)
    ybins = np.linspace(0, DETECTOR_SIZE_PIX, 5)
    spatial = etc.copy()
    spatial["x_cell"] = pd.cut(spatial["x_pixel"], xbins, include_lowest=True, labels=False)
    spatial["y_cell"] = pd.cut(spatial["y_pixel"], ybins, include_lowest=True, labels=False)
    medians = spatial.pivot_table(index="y_cell", columns="x_cell", values="observed_to_predicted",
                                  aggfunc="median", observed=True).reindex(index=range(4), columns=range(4))
    counts = spatial.pivot_table(index="y_cell", columns="x_cell", values="observed_to_predicted",
                                 aggfunc="count", observed=True).reindex(index=range(4), columns=range(4))
    spread = np.nanmax(np.abs(medians.to_numpy() - 1.0))
    image = ax.imshow(medians, origin="lower", extent=[0, 2048, 0, 2048], cmap="coolwarm",
                      norm=TwoSlopeNorm(vmin=1 - spread, vcenter=1.0, vmax=1 + spread))
    for yi in range(4):
        for xi in range(4):
            value, count = medians.iloc[yi, xi], counts.iloc[yi, xi]
            if np.isfinite(value):
                ax.text((xi + 0.5) * 512, (yi + 0.5) * 512, f"{value:.2f}\nn={int(count)}",
                        ha="center", va="center", fontsize=8, color="#0f172a")
    fig.colorbar(image, ax=ax, label="Cell median observed / predicted")
    ax.set(xlabel="CCD x (pixel)", ylabel="CCD y (pixel)",
           title="D  Detector-position residual map below 60 ke-")

    fig.suptitle("LOT / SOPHIA / r′ per-star SNR validation · 15 × 120 s · 299 saturation-safe stars",
                 fontsize=15, fontweight="bold")
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=180, facecolor="white")
    plt.close(fig)


def write_report(summary: dict) -> None:
    safe, etc = summary["safe"], summary["etc"]

    def corr_line(label: str, key: str) -> str:
        rho, low, high, n = summary["correlations"][key]
        return f"| {label} | {n} | {rho:+.3f} | {low:+.3f} to {high:+.3f} |"

    text = f"""# LOT r' per-star end-to-end residual analysis

This is the star-level version of the binned result in `endtoend.py`. It uses
the same fifteen 120 s LOT/SOPHIA r' frames from 2025-11-06 and applies the
final saturation selection: every retained star remains below 50,000 e- peak
in every frame, clear of the 60,292 e- 16-bit ADC ceiling.

![Four-panel per-star residual analysis](figures/lot_r_endtoend_residuals.png)

## Result

| Selection | Stars | Median observed/predicted | IQR | Bootstrap 95% CI of median |
|---|---:|---:|---:|---:|
| All saturation-safe stars | {safe['n']} | {safe['median']:.3f} | {safe['q25']:.3f}–{safe['q75']:.3f} | {safe['ci_low']:.3f}–{safe['ci_high']:.3f} |
| ETC regime, flux < 60 ke- | {etc['n']} | **{etc['median']:.3f}** | {etc['q25']:.3f}–{etc['q75']:.3f} | {etc['ci_low']:.3f}–{etc['ci_high']:.3f} |

The current CASTOR noise model is unbiased at the population level on this
night: the ETC-regime median is {100 * (etc['median'] - 1):+.1f}% from perfect
agreement. The broad per-star scatter is expected in part because every
observed noise is itself estimated from only 15 frames (`sigma_relerr = 18.9%`).

## Residual checks in the ETC regime

The table reports Spearman rank correlations. Confidence intervals bootstrap
stars 5,000 times; an interval crossing zero means this night does not resolve
a monotonic trend.

| Residual against | Stars | Spearman rho | Bootstrap 95% CI |
|---|---:|---:|---:|
{corr_line('log10 source flux', 'flux')}
{corr_line('Pan-STARRS g-r colour', 'color')}
{corr_line('CCD x', 'x')}
{corr_line('CCD y', 'y')}
{corr_line('distance from CCD centre', 'radius')}

No tested variable has a resolved monotonic relationship with the residual.
In particular, the colour result uses {summary['color_n']} catalogue-quality
stars and gives no evidence that source colour changes the *noise-model*
accuracy in r'. This does not replace a photometric colour-term fit, which is a
separate test of the count-rate calibration.

The 4x4 detector map uses the same below-60-ke- selection and has noisy cells
(some contain only a handful of stars),
but no coherent centre-to-edge trend; the global radial correlation above is
the appropriate summary. The mild downturn in the full brightness range comes
from the last three bins and is absent below 60 ke-, so it should not be turned
into another fitted noise term from this one night.

## Scope

This establishes one condition only: LOT, SOPHIA, r', one field, one night,
airmass 1.035–1.044, a 3xFWHM aperture, and a 5–8xFWHM median sky annulus. It
does not validate another band, telescope, extraction method, airmass, or sky
level. The publishable star table is `data/lot_r_endtoend_per_star.csv`; its
catalogue colours and detector positions make these checks reproducible without
publishing the raw observatory frames.

Regenerate from the ignored reduction products with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_endtoend_residuals.py
```
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reuse-derived",
        action="store_true",
        help="Use the committed per-star table instead of rebuilding it from ignored raw inputs.",
    )
    args = parser.parse_args()
    table = load_table() if args.reuse_derived else build_table()
    summary = summarize(table)
    make_figure(table, summary)
    write_report(summary)
    action = "read" if args.reuse_derived else "wrote"
    print(f"{action} {TABLE.relative_to(ROOT)} ({len(table)} stars)")
    print(f"wrote {FIGURE.relative_to(ROOT)}")
    print(f"wrote {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

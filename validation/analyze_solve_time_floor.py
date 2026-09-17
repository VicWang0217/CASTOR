"""Audit solve-for-time after the correlated background-flatness term was added.

``solve_required_exposures`` assumes SNR grows as sqrt(N). A flat-field or
background-gradient residual is correlated across a stack and instead creates
an asymptotic SNR ceiling. This script maps the discrepancy for the current
Lulin LOT and SLT presets without changing the engine.

Run from the repository root::

    uv run --with pandas --with matplotlib python validation/analyze_solve_time_floor.py
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from castor import physics, schema
from castor.calculator import run_calculation
from castorCLI import presets


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "validation/data/solve_time_flatness_audit.csv"
FIGURE = ROOT / "validation/figures/solve_time_flatness_floor.png"
REPORT = ROOT / "validation/report_sections/solve_time.md"

RIGS = {"LOT": "Sophia", "SLT": "SLT_DU934P"}
FILTERS = {"g": "Sloan_g", "r": "Sloan_r", "i": "Sloan_i"}
MAGNITUDES = np.arange(17.0, 23.01, 0.5)
TARGET_SNRS = (5.0, 10.0, 20.0, 30.0, 50.0)
SINGLE_EXPOSURE_S = 120.0
OBSERVING_TIME = "2026-01-15T16:00:00Z"
TARGET_RA_DEG = 113.65
TARGET_DEC_DEG = 31.89


def request_for(telescope: str, camera: str, filter_id: str, magnitude: float, target_snr: float):
    data = presets.load().resolve(
        "lulin", telescope=telescope, camera=camera, optic_filter=filter_id
    )
    data["instrument"]["throughput_correction"] = 1.0
    data["target"] = {
        "morphology": {"type": "point"},
        "brightness": {"type": "ab_mag", "target_mag": magnitude},
        "sed": {"type": "flat"},
        "ra": TARGET_RA_DEG,
        "dec": TARGET_DEC_DEG,
    }
    data["environment"].update({
        "observing_time_utc": OBSERVING_TIME,
        "auto_calc_background": False,
        "seeing_fwhm": 1.4,
        "diffraction_fwhm": 0.2,
        "optical_fwhm": 0.1,
        "tracking_fwhm": 0.1,
    })
    data["options"] = {
        "type": "solve_time",
        "aperture_factor": 0.85,
        "single_exp_time": SINGLE_EXPOSURE_S,
        "sky_annulus": {"inner_factor": 3.0, "outer_factor": 5.0, "estimator": "median"},
        "target_snr": target_snr,
    }
    return schema.ObservationRequest.model_validate(data)


def exact_exposure_count(request: schema.ObservationRequest, response, target_snr: float):
    """Solve the stack equation including its non-averaging flatness term."""
    instrument = request.instrument
    diagnostics = response.diagnostics
    budget = response.budget
    exposure = request.options.single_exp_time
    signal_frame = budget.source_count_rate * exposure
    background_pixels = diagnostics.num_pixels_aperture + diagnostics.num_pixels_sky_estimate
    random_variance_frame = (
        signal_frame
        + background_pixels
        * (
            budget.sky_count_rate * exposure
            + instrument.camera.dark_current_rate * exposure
            + instrument.camera.readout_noise**2
        )
    )
    flatness_amplitude_frame = (
        instrument.camera.background_flatness_fraction
        * budget.sky_count_rate
        * exposure
        * diagnostics.num_pixels_aperture
    )
    ceiling = math.inf if flatness_amplitude_frame == 0 else signal_frame / flatness_amplitude_frame

    denominator = signal_frame**2 - target_snr**2 * flatness_amplitude_frame**2
    if denominator <= 0:
        return ceiling, math.inf, math.nan
    exact = math.ceil(target_snr**2 * random_variance_frame / denominator)
    achieved = physics.calculate_total_snr(
        budget.source_count_rate,
        budget.sky_count_rate,
        instrument.camera.dark_current_rate,
        instrument.camera.readout_noise,
        diagnostics.num_pixels_aperture,
        exposure,
        exact * exposure,
        exact,
        diagnostics.num_pixels_sky_estimate,
        instrument.camera.background_flatness_fraction,
    )
    return ceiling, exact, float(achieved)


def build_table():
    rows = []
    for telescope, camera in RIGS.items():
        for band, filter_id in FILTERS.items():
            for magnitude in MAGNITUDES:
                for target_snr in TARGET_SNRS:
                    request = request_for(telescope, camera, filter_id, magnitude, target_snr)
                    response = run_calculation(request)
                    ceiling, exact, exact_achieved = exact_exposure_count(request, response, target_snr)
                    rows.append({
                        "telescope": telescope,
                        "camera": camera,
                        "band": band,
                        "ab_magnitude": magnitude,
                        "target_snr": target_snr,
                        "current_exposures": response.core.required_exposures,
                        "current_achieved_snr": response.core.total_snr,
                        "current_achieved_fraction": response.core.total_snr / target_snr,
                        "asymptotic_snr_ceiling": ceiling,
                        "exact_exposures": exact,
                        "exact_achieved_snr": exact_achieved,
                        "single_snr": response.core.single_snr,
                        "source_rate_e_s": response.budget.source_count_rate,
                        "sky_rate_e_s_pix": response.budget.sky_count_rate,
                        "aperture_pixels": response.diagnostics.num_pixels_aperture,
                        "sky_estimate_pixels": response.diagnostics.num_pixels_sky_estimate,
                        "background_flatness_fraction": request.instrument.camera.background_flatness_fraction,
                    })
    table = pd.DataFrame(rows)
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False, float_format="%.8g")
    return table


def make_figure(table):
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    colours = {"g": "#16a34a", "r": "#dc2626", "i": "#7c3aed"}
    slt = table[table["telescope"].eq("SLT")]

    ax = axes[0, 0]
    for band in FILTERS:
        subset = slt[(slt["band"].eq(band)) & (slt["target_snr"].eq(20))]
        ax.plot(subset["ab_magnitude"], subset["current_achieved_fraction"], "o-",
                color=colours[band], label=f"{band}'")
    lot = table[(table["telescope"].eq("LOT")) & (table["band"].eq("r"))
                & (table["target_snr"].eq(20))]
    ax.plot(lot["ab_magnitude"], lot["current_achieved_fraction"], "--",
            color="#0f172a", label="LOT r' control (f=0)")
    ax.axhline(1, color="#0f172a", linewidth=1, linestyle=":")
    ax.set(xlabel="Target AB magnitude", ylabel="Returned SNR / requested SNR",
           title="A  Current solve-for-time misses its own target")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ceiling = slt[slt["target_snr"].eq(5)]
    for band in FILTERS:
        subset = ceiling[ceiling["band"].eq(band)]
        ax.plot(subset["ab_magnitude"], subset["asymptotic_snr_ceiling"], "o-",
                color=colours[band], label=f"{band}'")
    for goal in (5, 10, 20):
        ax.axhline(goal, color="#64748b", linewidth=0.9, linestyle="--")
        ax.text(MAGNITUDES[-1] + 0.05, goal, f"SNR {goal}", va="center", fontsize=8)
    ax.set(yscale="log", xlabel="Target AB magnitude", ylabel="Asymptotic SNR ceiling",
           title="B  SLT's 2% correlated background term sets a ceiling")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    at_twenty = slt[slt["ab_magnitude"].eq(20)]
    for band in FILTERS:
        subset = at_twenty[at_twenty["band"].eq(band)]
        ax.plot(subset["target_snr"], subset["current_exposures"], "o-",
                color=colours[band], label=f"{band}' current")
        possible = subset[np.isfinite(subset["exact_exposures"])]
        ax.plot(possible["target_snr"], possible["exact_exposures"], "s--",
                color=colours[band], alpha=0.75, label=f"{band}' exact")
        impossible = subset[~np.isfinite(subset["exact_exposures"])]
        if len(impossible):
            ax.scatter(impossible["target_snr"], impossible["current_exposures"], marker="x",
                       s=75, linewidths=2, color=colours[band])
    ax.set(yscale="log", xlabel="Requested SNR for AB=20", ylabel="Number of 120 s exposures",
           title="C  Finite answers are returned above the ceiling")
    ax.legend(frameon=False, ncol=2, fontsize=8)

    ax = axes[1, 1]
    for band in FILTERS:
        row = at_twenty[(at_twenty["band"].eq(band)) & (at_twenty["target_snr"].eq(20))].iloc[0]
        request = request_for("SLT", "SLT_DU934P", FILTERS[band], 20, 20)
        response = run_calculation(request)
        n = np.unique(np.logspace(0, 5, 240).astype(int))
        snr = physics.calculate_total_snr(
            response.budget.source_count_rate,
            response.budget.sky_count_rate,
            request.instrument.camera.dark_current_rate,
            request.instrument.camera.readout_noise,
            response.diagnostics.num_pixels_aperture,
            SINGLE_EXPOSURE_S,
            n * SINGLE_EXPOSURE_S,
            n,
            response.diagnostics.num_pixels_sky_estimate,
            request.instrument.camera.background_flatness_fraction,
        )
        ax.plot(n * SINGLE_EXPOSURE_S / 3600, snr, color=colours[band], label=f"{band}'")
        ax.scatter(row["current_exposures"] * SINGLE_EXPOSURE_S / 3600,
                   row["current_achieved_snr"], color=colours[band], s=35)
    ax.axhline(20, color="#0f172a", linewidth=1, linestyle="--", label="requested SNR 20")
    ax.set(xscale="log", yscale="log", xlabel="Total integration (hours)", ylabel="Actual stacked SNR",
           title="D  More frames approach, but cannot cross, the floor")
    ax.legend(frameon=False)

    fig.suptitle("SLT solve-for-time audit under the measured 2% background-flatness floor",
                 fontsize=15, fontweight="bold")
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=180, facecolor="white")
    plt.close(fig)


def write_report(table):
    example = table[
        table["telescope"].eq("SLT")
        & table["ab_magnitude"].eq(20)
        & table["target_snr"].eq(20)
    ].set_index("band")

    rows = []
    for band in FILTERS:
        row = example.loc[band]
        exact = "unreachable" if not np.isfinite(row["exact_exposures"]) else str(int(row["exact_exposures"]))
        rows.append(
            f"| {band}' | {row['asymptotic_snr_ceiling']:.2f} | {int(row['current_exposures'])} | "
            f"{row['current_achieved_snr']:.2f} | {exact} |"
        )
    table_rows = "\n".join(rows)
    text = f"""# Solve-for-time versus the correlated background floor

The SLT DU934P preset now carries the 2% background-flatness residual measured
from real extended-source data. CASTOR correctly includes that non-averaging
term when it computes a stack's SNR, but still solves the number of exposures
with `N = (target_snr / single_snr)^2`. That square-root law is valid only when
every variance term is independent between frames.

![Solve-for-time flatness audit](figures/solve_time_flatness_floor.png)

## Concrete result

Standard case: SLT/DU934P, AB=20 point source, 120 s frames, 1.4 arcsec seeing,
0.85xFWHM aperture, 3–5xFWHM median annulus, the Lulin preset's own sky, target
near zenith, and requested SNR 20.

| Band | Model's asymptotic SNR ceiling | Current answer (frames) | SNR actually returned | Correct answer |
|---|---:|---:|---:|---:|
{table_rows}

The r' call is self-contradictory in one response: it says six exposures are
required and reports total SNR {example.loc['r', 'current_achieved_snr']:.2f},
below the requested 20. The i' target is farther beyond its ceiling. The g'
target is reachable, but needs {int(example.loc['g', 'exact_exposures'])}
frames rather than {int(example.loc['g', 'current_exposures'])}; the current
answer reaches only {example.loc['g', 'current_achieved_snr']:.2f}.

LOT is the control: Sophia's preset has `background_flatness_fraction = 0`, so
the old square-root law remains exact and its achieved/requested curve never
falls below one. Bright cases can overshoot because one indivisible frame
already exceeds the requested SNR.

## Correct algebra

For one frame, let `S` be source electrons, `V` the sum of every independent
variance term, and `F` the correlated flatness-noise amplitude. A stack of N
frames has

`SNR(N) = N*S / sqrt(N*V + N^2*F^2)`

and therefore the ceiling `S/F`. For a requested SNR `Q`:

`N = Q^2*V / (S^2 - Q^2*F^2)`

If the denominator is zero or negative, no finite exposure count can reach the
request under the model. The analysis table evaluates both the current and the
correct expression over LOT/SLT, g'/r'/i', AB 17–23 and target SNR 5–50.

## Consequence

This should be fixed before using solve-for-time with any camera whose
`background_flatness_fraction` is nonzero. The forward SNR calculation is
internally consistent; only the inverse solver assumes the superseded noise
law. A strict expected-failure test records the contradiction without silently
changing the engine as part of this analysis.

Regenerate with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_solve_time_floor.py
```
"""
    REPORT.write_text(text, encoding="utf-8")


def main():
    table = build_table()
    make_figure(table)
    write_report(table)
    print(f"wrote {TABLE.relative_to(ROOT)}, {FIGURE.relative_to(ROOT)}, and {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

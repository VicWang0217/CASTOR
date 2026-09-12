"""Regression for the LOT r' colour term measured from the public star table."""

import csv
import pathlib

import numpy as np


DATA = pathlib.Path(__file__).parent / "data" / "lot_r_endtoend_per_star.csv"


def _fit():
    with DATA.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = [
        row for row in rows
        if row["ps1_g_minus_r"]
        and float(row["snr_observed"]) >= 20
        and float(row["ps1_r_kron_mag"]) > 0
        and float(row["ps1_r_kron_mag_error"]) > 0
        and -0.20 <= float(row["ps1_r_mag"]) - float(row["ps1_r_kron_mag"]) <= 0.05
    ]
    colour = np.array([float(row["ps1_g_minus_r"]) for row in selected])
    zero_point = np.array([
        float(row["ps1_r_mag"])
        + 2.5 * np.log10(float(row["flux_e"]) / float(row["exposure_seconds"]))
        for row in selected
    ])
    reference = np.median(colour)
    keep = np.ones(len(colour), dtype=bool)
    for _ in range(10):
        design = np.column_stack([np.ones(keep.sum()), colour[keep] - reference])
        coefficients = np.linalg.lstsq(design, zero_point[keep], rcond=None)[0]
        residual = zero_point - (coefficients[0] + coefficients[1] * (colour - reference))
        centre = np.median(residual[keep])
        scatter = 1.4826 * np.median(np.abs(residual[keep] - centre))
        new_keep = np.abs(residual - centre) < 3 * scatter
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return selected, colour, zero_point, reference, coefficients, keep


def test_colour_term_uses_a_clean_point_source_sample():
    selected, _, _, _, _, keep = _fit()
    assert len(selected) == 101
    assert 90 <= int(keep.sum()) <= 100


def test_r_colour_term_is_small_but_positive():
    _, colour, zero_point, reference, coefficients, keep = _fit()
    residual = zero_point[keep] - (
        coefficients[0] + coefficients[1] * (colour[keep] - reference)
    )
    assert 0.0 < coefficients[1] < 0.05
    assert np.std(residual, ddof=2) < 0.02
    assert 23.75 < coefficients[0] < 23.82

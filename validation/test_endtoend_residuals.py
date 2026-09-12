"""Integrity checks for the publishable LOT r' per-star residual table."""

import csv
import pathlib

import numpy as np
import pytest

import endtoend
from castor import physics
from castorCLI import presets


DATA = pathlib.Path(__file__).parent / "data" / "lot_r_endtoend_per_star.csv"


@pytest.fixture(scope="module")
def stars():
    with DATA.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_table_is_the_final_saturation_safe_sample(stars):
    assert len(stars) == 299
    assert len({row["star_id"] for row in stars}) == len(stars)
    assert max(float(row["peak_max_e"]) for row in stars) < 50_000
    assert {int(float(row["n_frames"])) for row in stars} == {15}


def test_flux_bins_are_the_committed_end_to_end_measurement(stars):
    counts = {name: 0 for name in endtoend.MEASURED}
    for row in stars:
        counts[row["flux_bin"]] += 1
    assert counts == {name: values["n"] for name, values in endtoend.MEASURED.items()}


def test_stored_prediction_is_the_live_castor_noise_model(stars):
    camera = presets.load().profile("lulin").cameras["Sophia"].camera
    for row in stars:
        npix = float(row["aperture_pixels"])
        fwhm = float(row["fwhm_pixels"])
        exposure = float(row["exposure_seconds"])
        n_est = physics.calculate_sky_estimate_pixels(
            npix, 5.0, 8.0, fwhm, 1.0, estimator="median"
        )
        predicted = physics.calculate_single_snr(
            source_count_rate=float(row["flux_e"]) / exposure,
            sky_count_rate=float(row["sky_e_per_pixel"]) / exposure,
            dark_current_rate=camera.dark_current_rate,
            readout_noise=camera.readout_noise,
            num_pixels_aperture=npix,
            single_exp_time=exposure,
            num_pixels_sky_estimate=n_est,
            background_flatness_fraction=camera.background_flatness_fraction,
        )
        assert float(row["snr_predicted"]) == pytest.approx(predicted, rel=2e-6)


def test_etc_regime_is_unbiased_and_has_colour_and_position_coverage(stars):
    etc = [row for row in stars if row["in_etc_regime"] == "True"]
    ratios = np.array([float(row["observed_to_predicted"]) for row in etc])
    colours = [row["ps1_g_minus_r"] for row in etc if row["ps1_g_minus_r"]]
    assert len(etc) == 261
    assert len(colours) >= 250
    assert np.median(ratios) == pytest.approx(1.0, abs=0.02)
    assert all(0 <= float(row["x_pixel"]) < 2048 for row in stars)
    assert all(0 <= float(row["y_pixel"]) < 2048 for row in stars)

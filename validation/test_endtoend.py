"""Does CASTOR's own noise model reproduce the scatter `endtoend.py` measured?

The prediction is recomputed here from the shipped Sophia preset and
`endtoend.py`'s committed conditions, rather than stored as a fixed number --
a change to `physics.py`, or to presets.json's read noise or dark current,
changes what this test expects, not just what it observes.
"""
import pytest

import endtoend
from castor import physics
from castorCLI import presets


@pytest.fixture(scope="module")
def camera():
    return presets.load().profile("lulin").cameras["Sophia"].camera


def _predicted_snr(camera, flux_e):
    """CASTOR's own single-exposure SNR for `flux_e` electrons under
    `endtoend.CONDITIONS`, including the sky-estimate term the 2026-08-30
    check found missing (ATBD 4.3.2 / 5.3, "Cost of the Sky Estimate")."""
    c = endtoend.CONDITIONS
    n_est = physics.calculate_sky_estimate_pixels(
        num_pixels_aperture=c["num_pixels_aperture"],
        inner_factor=5.0, outer_factor=8.0,
        # Already in pixels, so pixel_scale=1 leaves the geometry untouched --
        # the annulus/aperture ratio is all that matters here, not arcsec.
        total_fwhm=c["fwhm_pix"], pixel_scale=1.0,
    )
    return physics.calculate_single_snr(
        source_count_rate=flux_e / c["exptime_s"],
        sky_count_rate=c["sky_e_pix"] / c["exptime_s"],
        dark_current_rate=camera.dark_current_rate,
        readout_noise=camera.readout_noise,
        num_pixels_aperture=c["num_pixels_aperture"],
        single_exp_time=c["exptime_s"],
        num_pixels_sky_estimate=n_est,
    )


#: How far obs/pred may sit from 1.0 in each bin. Two things widen it: a
#: single representative point standing in for the whole bin's flux range
#: (worst at the faint end, where the SNR-flux curve bends fastest), and
#: `sigma_relerr` -- the uncertainty on a scatter measured from 15 frames,
#: 19% for any one star -- shrinking only as sqrt(n) with bin size.
TOLERANCE = {
    "<3k": 0.10, "3-10k": 0.08, "10-30k": 0.08, "30-60k": 0.05,
    "60-120k": 0.08, "120-240k": 0.10, ">240k": 0.10,
}


@pytest.mark.parametrize("flux_bin", list(endtoend.MEASURED))
def test_predicted_snr_matches_measured_scatter(camera, flux_bin):
    """If this breaks, the noise model moved -- not this measurement.

    `endtoend.MEASURED` is fixed, from one real night; the prediction is
    recomputed from the live preset and physics module every run.
    """
    m = endtoend.MEASURED[flux_bin]
    ratio = m["snr_obs"] / _predicted_snr(camera, m["flux_e"])
    assert ratio == pytest.approx(1.0, abs=TOLERANCE[flux_bin])


def test_the_saturation_cut_still_sits_below_the_real_ceiling(camera):
    """`endtoend.py`'s 50000 e- peak cut must stay below the camera's real
    saturation -- 65535 ADU x 0.92 e-/ADU for this gain setting -- not
    presets.json's `full_well_capacity`, which the delivered, gain-corrected
    frames never get near. See `data/raw/_endtoend_2026-08-30/RESULT.md`'s
    2026-09-10 follow-up. If SOPHIA's delivered bit depth or gain ever
    changes, this is the number to revisit before trusting the >240k bin
    again.
    """
    adc_ceiling_e = 65535 * 0.92
    assert adc_ceiling_e > 50_000
    assert camera.full_well_capacity > adc_ceiling_e


@pytest.mark.parametrize("flux_bin", list(endtoend.EXTRACTION_STABILITY))
def test_tight_aperture_excess_is_extraction_not_intrinsic(flux_bin):
    """The 0.85 x FWHM shortfall must not become a source-noise term.

    These are the same stars on the same frames. Intrinsic variability would
    survive both the free-width PSF fit and the wider aperture; instead those
    two agree while only the tight aperture has extra scatter. The 0.85
    default is consequently a photon-limit operating point, not a promise that
    uncorrected aperture photometry will attain that limit.
    """
    measured = endtoend.EXTRACTION_STABILITY[flux_bin]
    baseline = max(measured["aperture_150_rms"], measured["psf_rms"])
    assert measured["aperture_085_rms"] > 1.5 * baseline
    assert measured["aperture_150_rms"] == pytest.approx(
        measured["psf_rms"], abs=0.003
    )

"""LOT/SOPHIA against real frame-to-frame scatter: does CASTOR predict the
noise a night of imaging actually delivers?

Every other file under validation/ checks one of CASTOR's inputs -- throughput,
sky brightness, extinction. This is the only one that closes the loop: it
compares CASTOR's *predicted SNR* against *measured* scatter, on data where
the signal is not in question -- LOT's throughput was fitted from these same
frames, so the electron counts agree by construction, and only the noise model
is on trial.

Method, so the numbers can be argued with:

* LOT/SOPHIA, r', night of 2025-11-06: 15 frames of 120 s, airmass 1.035-1.044,
  FWHM 2.3-2.6 pix. 314 Pan-STARRS stars measured on every frame, 3xFWHM
  aperture, 5-8xFWHM sky annulus (median), isolated by 8xFWHM.
* `snr_obs` is the star's own scatter across the fifteen frames -- not a
  model -- after dividing out the 0.30% rms transparency drift common to every
  star in a frame.
* Stars whose peak pixel exceeds 50000 e- in any frame are excluded. The
  camera's real ceiling at this gain is 65535 ADU x 0.92 e-/ADU = 60292.2 e- --
  the 16-bit ADC, not the 150000 e- `full_well_capacity` in presets.json, which
  this delivered, gain-corrected data never gets near. 15 of the original 314
  stars clip at exactly that number or ride the non-linear roll-off just below
  it; without the cut the brightest bin's obs/pred collapses to 0.22 for a
  reason that has nothing to do with the noise model. See
  `data/raw/_endtoend_2026-08-30/RESULT.md`'s 2026-09-10 follow-up for how
  that was found, and `test_endtoend.py` for the ceiling this relies on.

The frames (~2 GB) and the scripts that reduced them stay out of the
repository, in `data/raw/_endtoend_2026-08-30/` (gitignored). What is
committed is the result of that reduction, below -- run `test_endtoend.py`
without needing any of it.
"""

#: Aperture geometry and sky conditions, shared across every bin below: one
#: field, one night, one aperture, so these move less than 0.3% bin to bin.
#: Feeds `test_endtoend.py`'s recomputation of CASTOR's own prediction.
CONDITIONS = dict(
    sky_e_pix=1577.9,           # e-/pix over the exposure, from the annulus
    num_pixels_aperture=169.8,
    fwhm_pix=2.4494,            # feeds calculate_sky_estimate_pixels
    exptime_s=120.0,
)

#: Measured 2026-08-30, corrected 2026-09-10 for the saturation cut above.
#: One row per flux bin: how many of the 314 stars land in it, the median
#: source signal, and the median measured SNR. The bottom four bins are
#: `<60k`'s narrower sub-ranges from the original write-up rather than one
#: combined bin -- a single point can't represent the SNR-vs-flux curve
#: across a 20x span, so a wide "<60k" bin does not reproduce the per-star
#: median it is meant to summarise (it read 1.14, not the ~1.0 every
#: sub-range actually lands on).
#:
#:   n         stars in the bin, after the 50000 e- peak cut
#:   flux_e    median source electrons in the aperture over 120 s
#:   snr_obs   median measured SNR: flux_e / (scatter across the 15 frames)
MEASURED = {
    "<3k":      dict(n=69,  flux_e=2162.7,   snr_obs=3.26),
    "3-10k":    dict(n=112, flux_e=4517.3,   snr_obs=7.72),
    "10-30k":   dict(n=56,  flux_e=16853.7,  snr_obs=28.15),
    "30-60k":   dict(n=24,  flux_e=40151.4,  snr_obs=62.01),
    "60-120k":  dict(n=18,  flux_e=84235.6,  snr_obs=127.57),
    "120-240k": dict(n=15,  flux_e=177261.0, snr_obs=219.00),
    ">240k":    dict(n=5,   flux_e=261854.0, snr_obs=342.48),
}

#: HAP-11 follow-up on the same fifteen frames. Each entry is the median
#: per-star fractional RMS after removing the transparency change common to a
#: frame. The comparison uses only stars present on every frame, below the
#: measured 50000 e- linearity cut, and with an unflagged elliptical-Gaussian
#: PSF fit on every frame. Flux bins use the 3 x FWHM aperture signal.
#:
#: A real intrinsic change in a star is present whichever extractor measures
#: it. Instead, the free-width PSF fit and the 1.5 x FWHM aperture agree while
#: the same stars measured at 0.85 x FWHM carry substantially more scatter.
#: That isolates the excess as tight-aperture extraction sensitivity, not a
#: source noise term for CASTOR to add. The ignored raw-data reduction is
#: `data/raw/_endtoend_2026-08-30/psf_scan.py`.
EXTRACTION_STABILITY = {
    "30-60k":  dict(n=24, aperture_085_rms=0.033286, aperture_150_rms=0.009289,
                    psf_rms=0.009702),
    "60-120k": dict(n=17, aperture_085_rms=0.012708, aperture_150_rms=0.006452,
                    psf_rms=0.007170),
    ">120k":   dict(n=20, aperture_085_rms=0.015579, aperture_150_rms=0.004076,
                    psf_rms=0.005845),
}

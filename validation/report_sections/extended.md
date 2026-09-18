# Extended-source noise: model against a real galaxy

The extended-source noise model was checked against NGC 3621 in SLT r'
(42 x 30 s, 2024-04-12, Siril mean-stack). The measured noise carries no model
at all: the clear frames are split alternately into two half-stacks and
differenced, so nothing astrophysical survives and the aperture scatter of the
difference is noise and noise only. This is the committed, publishable summary
of the fuller reduction whose raw frames stay out of the repository.

![Extended-source noise vs aperture radius](figures/extended_noise_vs_radius.png)

## What the picture shows

The left panel is aperture noise against aperture radius, scaled to a single
30 s frame, for six radii from 1.5" to 12.2". The measured curve (black) is
compared with what CASTOR predicted at the time of the test and what it
predicts after the fix on branch `fix/slt-noise-model`:

* **Shipped model (red), Andor read noise 3.3 e-, no correlated term.** It
  grows only as sqrt(area) and falls progressively below the measurement, from
  1.5x low at 1.5" to 3.2x low at 12.2".
* **Fixed model (blue), read noise 9.28 e- plus a 2% background-flatness
  floor.** It tracks the measurement to about 1% once the flatness term
  dominates past ~4", and stays within 6-10% at the smallest apertures where
  the read-noise correction alone carries the change.

The right panel is the same story as a ratio: the shipped prediction runs away
from unity with radius, the fixed prediction sits on it.

## Two independent fixes are folded into the blue curve

**Read noise.** The delivered calibrated frames carry ~9.3 e-/px of non-Poisson
noise, measured twice independently, not the 3.3 e- the preset held for the
DU934P. That is the noise a user of these pipeline products actually gets.

**Background flatness.** All six radii fit `sigma^2 = N_pix*(sky + RN^2) +
(f*N_pix)^2` with `f = 0.741 e-`, 2.0% of the background. It enters as area,
not sqrt(area) -- flat-field and background-gradient residuals rather than
photons -- so past ~4" it, not read noise, is what the gap was made of. The
floor measured on this cloudy night is an upper limit on a good night.

## Scope

One telescope, one band, one galaxy, one night. The rate formula itself is
correct to 1 part in 10^4; this tests only the noise and the aperture geometry.
Two known design limitations are recorded separately and not visible here: the
aperture is tied to the PSF rather than specifiable in arcsec, and
`saturation_time_limit` assumes a flat profile.

Regenerate the figure from the committed per-radius aggregate with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_extended_noise.py
```

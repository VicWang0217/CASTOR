# Summary

Five analyses now connect the real Lulin frames to decisions an observer can
make with CASTOR. Each result has a publishable table, a regenerating script,
and a validation test; the observatory's FITS files remain ignored.

| Analysis | Main result | Section |
|---|---|---|
| LOT r' per-star SNR | 261 stars below 60 ke- give median observed/predicted SNR **1.011**, bootstrap 95% CI 0.980–1.050; no resolved flux, colour, or detector-radius trend | End-to-end residuals |
| LOT r' colour term | `+0.0182 mag per mag` in Pan-STARRS g-r (95% CI +0.0036 to +0.0423); only 1.5% across the central 90% colour range, far too small to explain the r' throughput excess | Colour term |
| Correlated-noise inverse audit | SLT r', AB=20, target SNR 20 returns six frames but only SNR **14.44**; the model ceiling is 18.46, so the request is unreachable | Solve-for-time flatness |
| Forward performance | Under the stated standard scene, one-hour SNR=5 limits are LOT g/r/i = **23.91/23.84/23.04**, SLT = **21.84/21.36/20.51** AB mag | Forward performance |
| Extended-source noise | NGC 3621 SLT r': the fixed model (RN 9.28 e- + 2% flatness) tracks measured aperture noise to **~1%** to 12"; the shipped RN-3.3 model ran 1.5–3.2x low | Extended-source noise |

Four further figures draw data the tables only summarise: the
observed-vs-predicted noise 1:1 scatter and the ADC-ceiling saturation plot
(`figures/lot_r_noise_scatter.png`), the per-term noise budget with depth-vs-
time and the validated anchor (`figures/lulin_noise_budget.png`), the
solve-for-time reachability map (`figures/solve_time_reachability.png`), and
the extended-source noise-vs-radius comparison
(`figures/extended_noise_vs_radius.png`).

## What is established

The current forward noise equation reproduces the LOT r' frame-to-frame
scatter at the population level after the real ADC saturation selection is
applied. Within the statistical power of fifteen frames, its residual does not
depend monotonically on source flux, stellar colour, CCD x/y, or distance from
the detector centre.

LOT r' does have a small natural-system colour coefficient against Pan-STARRS,
but the reference-colour zero point agrees with the existing multi-night
calibration within 0.0033 mag. The present throughput therefore need not move
on this evidence. Future calibration should nevertheless reject extended
sources and fit colour explicitly.

## What is broken

The forward and inverse calculations disagree whenever
`background_flatness_fraction` is nonzero. Forward stacking correctly retains
the correlated term; `solve_required_exposures()` still assumes every term
averages down as sqrt(N). Question 17 and a strict xfail record the defect. The
forward-performance atlas avoids it by evaluating actual integer stacks.

## What the capability comparison does not prove

SLT carries a measured 2% correlated background residual. Sophia currently
carries zero, meaning unmeasured rather than absent, so the long-integration
LOT/SLT gap is partly a model-completeness gap. The atlas is also conditional
on one pointing, 1.4 arcsec seeing, 120 s frames, the top-hat passbands and the
site-wide extinction fallback. It is a planning baseline, not an empirical
limiting-magnitude claim.

## Recommended order from here

1. Fix the inverse correlated-noise algebra and represent unreachable SNR goals
   in the response schema, CLI and GUI.
2. Measure Sophia's background-flatness floor with the same empty-aperture
   method used for SLT before treating long-stack LOT/SLT curves as symmetric.
3. Repeat the colour-term fit in g' and i' and on independent fields.
4. Repeat the end-to-end SNR check at another sky level, band and airmass.
5. Add uncertainty propagation only after the two instruments' missing floors
   and colour terms are represented; otherwise the interval would formalise
   known omissions as if they were zero.

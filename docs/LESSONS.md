# Lessons

Traps this project has already fallen into, so it does not fall into them
twice. Each entry is a thing that once looked right and was not: what it looked
like, why it was actually wrong, how it surfaced, how to see it for yourself, and
what now stops it from coming back — with a commit or file you can follow.

This is not a to-do list. What is still *open* — the things nobody has decided
or measured yet — lives in [`../validation/QUESTIONS.md`](../validation/QUESTIONS.md),
the single place for that, with the same field on every row: who can close it.
A lesson lands here once it is closed and understood; if a row here still has an
open question hanging off it, that question belongs in QUESTIONS.md too, and this
entry links to it rather than restating it.

Each entry carries four lines. **Found by** — the check or moment that first
exposed it. **Reproduce** — how to see the evidence today; where the bug is
fixed, this is the regression test that now pins it, so a green run *is* the
proof it stays fixed. **Guard** — what stops it recurring. Commit hashes in
parentheses are the source.

---

## Physics that the units cannot catch

### The same erg/s/cm²/Å can be in two different frames of reference

The heaviest bug this engine has carried. All three count rates — star, galaxy,
sky — shared `_calculate_base_electron_rate`, so the sky was routed through the
same `10^(-0.4·k·X)` extinction term as target starlight. But `mu_sky` descends
from `mu_dark`, a brightness *measured from the ground*, through the very
atmosphere that coefficient describes. Attenuating it again counts the
atmosphere twice — and worse, inverts the sign: sky surface brightness *rises*
with airmass, and the engine was making it fall. Both sides were plain floats in
erg/s/cm²/Å, so the dimensions were right and the frame of reference was wrong
(`b1df2a7`).

- **Found by** building the validation suite against ESO's FORS2 ETC and LCO's
  calculator. FORS2 puts the sky 17.6% brighter at X=1.5; we were 6% fainter —
  nobody agreed with us in either direction.
- **Reproduce** `pytest validation/test_eso.py`. The airmass-growth case is a
  strict xfail (flat sky is only exact at the zenith); the sign and
  double-count regression is asserted alongside it.
- **Guard** a top-of-atmosphere flux (star, galaxy) and a ground-measured
  surface brightness (sky) are different quantities even sharing a unit. They now
  take different paths: `calculate_sky_background_rate` does not accept `airmass`
  or `extinction_coeff` at all — "a parameter taken and discarded is an
  invitation to pass it and assume it did something, which is how this got here."

### An unobservable target is reported as a very faint observable one

`calculator.py` clamps the zenith angle to 89°. A target 35° *below* the horizon
comes back not as an error but as airmass 57 — a real-looking number for a source
that cannot be seen (its true airmass is −1.8). The clamp keeps the maths finite
at the cost of making invalidity look like faintness (`44ab4af`).

- **Found by** pinning the CLI test times: the saturation test's target, twelve
  hours off transit, was not low but below the horizon, and the engine answered
  anyway.
- **Reproduce** ask the engine for a target well below the horizon at a fixed
  time; the returned airmass is a large positive number, not a rejection.
- **Guard** reachability is not brightness. Do not read a large airmass as "hard
  to observe" — check the target is actually up.

## Noise and calibration

### A missing noise term always looks like a calibration offset first

The noise model kept coming out optimistic against real frames, and the first
instinct each time was a throughput or zero-point tweak — one multiplicative
offset. It was never that. Twice the cause was a term the model did not charge
for: the sky estimate is not free (subtracting it adds variance scaling as
`k_ap²`, `051994e`), and the background is not perfectly flat (the NGC 3621
extended-source check was short a flatness term, `6a34861`).

- **Found by** the end-to-end SNR check (`validation/endtoend.py`) comparing
  predicted SNR to real frame-to-frame scatter.
- **Reproduce** `pytest validation`, then read `validation/ENDTOEND_RESIDUALS.md`
  — residuals resolved by flux.
- **Guard** when observed SNR sits below predicted, rule out a *missing noise
  source* before touching throughput or the zero point. An offset moves every
  flux bin by the same fraction; a missing noise term does not — that is how the
  two are told apart.

### The datasheet describes the detector; the observer measures the instrument

Two of the same shape. SOPHIA's read noise was kept at the datasheet's 1 MHz
7.0 e⁻ over a measured 7.9, on the sound reasoning that a PTC intercept bounds
the *sensor* — but an ETC predicts what the *observer* measures, and the observer
gets the residual pattern noise too (`2128f7d`). Separately, SLT's read noise was
a real DU934P datasheet number, 3.3 e⁻, for the wrong readout port (`f393a54`).
Both plausible, both passed review, both understated noise in the unsafe
direction.

- **Found by** a photon-transfer curve over 89 adjacent frame pairs (SOPHIA), and
  tracing SLT's value back to its datasheet row (which port).
- **Reproduce** compare the `read_noise` in `presets.json` against the
  frame-measured value recorded in its provenance row.
- **Guard** use the value the frames show, not the datasheet's, and record the
  readout configuration (port, speed, temperature) beside the number — every
  camera in the 2011 prototype spans ×2–12 across its ports.

### The reference data has limits too — check them before blaming the model

Above ~60 ke⁻ the observed SNR fell short and the shortfall grew with brightness
— the shape of a real model defect (closed HAP-72). It was not one: the reference
frames' own 16-bit ADC was saturating, so the "systematic" was in the data
CASTOR was judged against, not in CASTOR (`a4dd05c`).

- **Found by** binning the end-to-end residual by electrons in the aperture and
  watching it turn over only at the bright end.
- **Reproduce** `pytest validation`; the per-flux-bin table in
  `ENDTOEND_RESIDUALS.md` shows the turnover above the linear range.
- **Guard** a discrepancy that appears only at the bright end is a saturation
  suspect first. Establish the reference instrument's linear range before adding
  a term to chase a residual outside it.

### You cannot fit extinction from a night that was not photometric

Per-band extinction coefficients were fitted, shipped, and marked CLOSED — then
retracted. The night (2024-04-14) had an hour of cloud and never fully
recovered, and because the target was setting, transparency was fading while
airmass rose. `zp = zp0 − k·X` cannot separate those two moving together: the fit
charges the fade to airmass and inflates `k`, band by band in step with how much
each faded — which is what identifies it as weather, not a dusty site (`ad0c78b`).

- **Found by** a model-free test: frames at the *same* airmass early and late in
  the night have an identical extinction term by definition, so any zero-point
  gap between them is transparency alone — and the sky was fainter later in every
  band.
- **Reproduce** `pytest validation/test_slt.py`; it asserts the per-band values
  are ABSENT and pins the fade/excess correlation. `slt.WHY_NO_EXTINCTION` holds
  the numbers.
- **Guard** extinction needs a night *verified* photometric — returning to the
  same airmass and finding the same zero point, not merely spanning a range.

### A borrowed model needs one in-situ check before you trust it

Standing SkyCalc's Paranal sky in for Lulin looked reasonable until checked
against our own frames: Lulin is brighter in every band, and *monotonically
worst in the blue*. Airglow lives in the near-infrared, so an airglow excess
would be red-weighted — this is the other way round, which is scattered
artificial light and aerosol, a term SkyCalc has for neither (`5b45ad5`).

- **Found by** comparing SkyCalc's absolute moonless sky at our sightline against
  the absolute sky the 123 frames measure.
- **Reproduce** `pytest validation/test_sky_model.py` (and `validation/skycalc.py`
  for the component split).
- **Guard** when you borrow another site's model, check the one quantity you can
  against local data first. Only what is genuinely site-independent transfers —
  here, the interplanetary zodiacal component, anchored throughout to our own
  photometry.

## Software structure

### The same physics implemented twice will drift, silently

A band's throughput is keyed by telescope (a value measured on LOT says nothing
about SLT). The browser applied that fragment flat, naming a field
(`instrument.telescope.LOT`) that does not exist, and `applyFragment` silently
skips fields it cannot find — so the measured override applied as *nothing*, and
the rig's default stood in. The shipped page was 33% low in r' with nothing to
announce it, because the overrides beside it worked (`2d93107`).

- **Found by** the engine having a test for this exact bug (the schema fix) while
  the form had none — "which is why one was fixed and the other went on being
  wrong."
- **Reproduce** check the fragment's field path against the schema: a name the
  model does not define is dropped without error.
- **Guard** when the same physics lives in two places (engine and client), a
  divergence is invisible until a test pins them to the same answer. Give the
  second copy the first copy's test.

### An inverse solver must invert the same noise model the forward path uses

`calculate_total_snr()` correctly models the 2% background-flatness term as
correlated across a stack, creating an asymptotic SNR *ceiling*. But
`solve_required_exposures()` still returned `(target_snr / single_snr)²`, which
assumes every noise term averages down as √N. So the solver promised six frames
reach SNR 20 when its own forward model puts the ceiling below 20 and no frame
count gets there (`14e0e6c`).

- **Found by** the forward model and the inverse solver disagreeing on one case:
  SLT/DU934P, r', AB=20, 120 s frames, requested SNR 20 — six frames "enough",
  achieved SNR 14.44.
- **Reproduce** `pytest validation/test_solve_time_floor.py` (strict xfail pinning
  the contradiction); QUESTIONS.md item 17 records the fix.
- **Guard** a forward model and its inverse must share one noise model. A
  correlated term breaks the √N shortcut; the solver must invert
  `SNR(N)=N·S/√(N·V+N²·F²)` and return an explicit unreachable result past S/F.

### Prose and diagrams drift from the schema they describe

`architecture.md` claimed `environment` was "a discriminated union from the
start". It is not — `EnvironmentCondition` is a flat model with no discriminator;
the pillar that genuinely is a union is `options`, on `type`. The words and the
picture had drifted from the models and disagreed with each other (`090c40e`,
`f0ddcf7`).

- **Found by** reading the prose next to `schema.py` and finding they disagreed.
- **Reproduce** grep `schema.py` for a discriminator on `EnvironmentCondition` —
  there is none; the `Union`/discriminator is on the options model.
- **Guard** documentation is only true the moment it is written. When a claim
  names a structural property (a union, a required field, a count), check it
  against the code, not the last version of the prose. This file included.

## Packaging and robustness

### A planning tool must not crash on planning ahead

astropy's `iers.conf.auto_max_age` defaults to 30 days: once the bundled
Earth-orientation table is that stale relative to *now*, an AltAz transform does
not extrapolate, it raises — and a packaged desktop app has no network fallback
to fetch a fresh table. So any observation more than ~a month out raised
`ValueError`, in a calculator whose whole job is planning ahead (`dc15bbd`).

- **Found by** offline planning for a date more than a month out crashing in the
  desktop build.
- **Reproduce** call `moon.get_moon_and_target_geometry` for a date >30 days out
  with `astropy` downloads forced off; without the fix it raises. Confirmed
  against 2026-12-01 and 2028-06-01.
- **Guard** `auto_max_age = None` is set deliberately. It costs no real precision
  — astropy's out-of-coverage fallback is a 50-year polar-motion mean, good to the
  arcsecond, well inside any seeing FWHM this tool sees.

### A test that reads the clock passes until it doesn't

Every calculating test in `tests/test_cli.py` let the CLI default its observing
time to *now*, so each silently asked where the sky happened to be. They passed
when written and failed the same afternoon with nothing changed: the target had
set, and a source that filled the well in half a second was six thousand times
fainter twelve hours later (`44ab4af`).

- **Found by** the suite going red in the afternoon with no code change since the
  morning it was written.
- **Reproduce** historically, running `tests/test_cli.py` at different times of
  day; now every ephemeris-touching test passes an explicit `--time`, so it does
  not recur.
- **Guard** any test that touches the ephemeris pins a `--time`. A test whose
  result depends on wall-clock time is not testing what it claims to.

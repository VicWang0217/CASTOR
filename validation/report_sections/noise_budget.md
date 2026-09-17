# Noise budget and depth versus time

The forward atlas above, read for what limits each rig rather than for headline
numbers.

![Noise budget and depth vs time](figures/lulin_noise_budget.png)

The top panels split the one-hour variance by term against source magnitude.
LOT r' walks from source-limited at the bright end through a sky-dominated
middle to a read-noise floor; SLT r' is swamped by the correlated 2% flatness
term at everything fainter than roughly AB 17, which is exactly why its long
stacks stop getting deeper. The bottom panel is the consequence: SNR=5 limiting
magnitude against total exposure time, where LOT keeps descending as sqrt(t)
while the SLT curves bend flat. The star marks the 1,800 s LOT r' point where
the underlying noise model was validated end-to-end at obs/pred 1.011 (the
per-star end-to-end section), so this one abscissa is anchored to data and the
rest is model extrapolation.

Regenerate with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_noise_budget.py
```

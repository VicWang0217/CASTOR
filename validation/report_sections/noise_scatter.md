# Per-star noise scatter and the ADC ceiling

Two per-star views the binned end-to-end table cannot show, both read from the
committed `data/lot_r_endtoend_per_star.csv` and
`data/lot_r_endtoend_peak_vs_flux.csv`.

![Noise 1:1 scatter and the ADC ceiling](figures/lot_r_noise_scatter.png)

The left panel is the model-free validation: predicted against observed noise
for every star, with the 1:1 line. The faint majority pile up at the sky-plus-
read floor near 600 e-, scattering symmetrically about the line, and the
ETC-regime median lands at obs/pred 1.011. The right panel is the bright-end
story in one plot: each star's brightest pixel against its aperture flux,
flattening hard at the 60,292 e- 16-bit ADC ceiling. Six stars are pinned there
and another handful clear the 50,000 e- peak cut, which is why the headline
sample is restricted to saturation-safe stars rather than carrying a spurious
noise floor.

Regenerate with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_noise_scatter.py
```

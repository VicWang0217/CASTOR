# Solve-for-time reachability map

The single operating point in the section above is one cell of a larger grid.
Sweeping LOT/SLT x g'/r'/i' x AB 17-23 x requested SNR 5-50 shows where the
defect lives.

![Solve-for-time reachability map](figures/solve_time_reachability.png)

Colour is achieved-over-requested SNR, so 1.0 is honest, blue is a harmless
integer-frame overshoot, and red is a request the solver accepts but
undershoots. The LOT row (flatness 0) is honest everywhere because the old
sqrt(N) law is still exact there. The SLT row turns red across a wide band of
faint targets and high SNR goals, and the hatched cells are requests above the
model's asymptotic ceiling, unreachable at any exposure. Worst case is SLT i'
at AB 23, SNR 50, where the current answer delivers 1% of the request.

Regenerate with:

```bash
uv run --with pandas --with matplotlib python validation/analyze_solve_reachability.py
```

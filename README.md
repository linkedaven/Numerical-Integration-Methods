# Numerical Integration Methods Explorer

An interactive tool comparing four numerical integration methods at
arbitrary precision, with a live convergence plot and results table.

## Methods

1. **Trapezoidal rule** (composite, degree-1 Newton-Cotes)
2. **Simpson's 1/3 rule** (composite, degree-2 Newton-Cotes)
3. **Simpson's 3/8 rule** (composite, degree-3 Newton-Cotes)
4. **Gaussian quadrature** (Gauss-Legendre; nodes/weights computed at
   working precision via the Golub-Welsch eigenvalue method, not looked up
   from a fixed-precision table)

## Why arbitrary precision?

Gaussian quadrature converges so fast on a smooth function that with
ordinary double-precision (float64) numbers, its error hits the ~1e-16
machine-epsilon floor almost immediately. Past that point you're no longer
seeing the method's real truncation error -- just double-precision
round-off noise, which doesn't shrink monotonically and can even tick back
*up* slightly as more terms are summed.

This project uses [mpmath](https://mpmath.org/) at 60 decimal digits of
working precision, pushing that round-off floor down to roughly 1e-60 --
far below anything reached in this demo -- so every convergence curve shows
genuine truncation error across its full range, not an artifact of
floating-point precision.

## Features

- Interactive **n slider** (1-60) in the plot window itself, with a **Run**
  button to recompute
- Live convergence plot (log-log error vs. n) for all four methods
- Results table showing the approximation and error at the chosen n
- Dark theme UI, consistent with this author's other simulation projects
- `gauss_legendre_nodes_weights` and other expensive per-n computations are
  memoized with `functools.lru_cache`, so repeated use of the slider/Run
  button doesn't redo work already computed

## Requirements

- Python 3.9+
- See `requirements.txt`

## Setup

```bash
pip install -r requirements.txt
```

## Usage

```bash
python integral_sol.py
```

Drag the **n** slider to the number of points/subintervals you want, then
click **Run**. The convergence plot redraws showing every n from 1 up to
your chosen value, and the results table updates to show each method's
approximation and error at that n.

## Test problem

By default, the script evaluates:

```
I = integral of sin(x) from 0 to pi = 2
```

To try a different function, edit `f(x)` near the top of the script (it's
an `mpmath`-based function, so use `mp.sin`, `mp.exp`, etc. rather than
plain Python `math`), and update `a`, `b`, and `I_exact` to match your new
problem and its known closed-form answer.

## License

MIT -- see [LICENSE](LICENSE).

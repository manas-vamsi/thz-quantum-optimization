"""Model 2 -- Fibonacci-based antenna layouts.

IMPORTANT DISTINCTION (see docs/fibonacci_vs_golden.md)
-------------------------------------------------------
The golden-ratio logarithmic spiral uses the *irrational constant*
``PHI = (1+sqrt(5))/2`` inside a continuous radial law.  The Fibonacci sequence
``F_n = F_{n-1} + F_{n-2}`` is a *discrete integer sequence* whose successive
ratios converge to PHI.  Substituting ``F_n`` for ``PHI`` in the spiral equation
is not a defined operation, so this module implements three explicit,
mathematically stated Fibonacci constructions instead.  Each is selected with
``variant=``:

``"radial"``
    Fibonacci numbers set the *radii*:

        r_n = r_min + (r_max - r_min) * (F_n - F_min) / (F_max - F_min)
        theta_n = theta_0 + n * theta_step

    ``theta_step`` defaults to the golden angle but is fully configurable.
    Because ``F_n ~ PHI**n / sqrt(5)``, the normalised radial sequence is close
    to (but not identical to) a geometric progression; the departure is largest
    for the first few terms, where the Fibonacci sequence is not yet in its
    asymptotic regime.

``"golden_angle"``
    Vogel / phyllotaxis construction, the classical "Fibonacci sunflower":

        r_n = sqrt(n / N) scaled onto [r_min, r_max]
        theta_n = theta_0 + n * GOLDEN_ANGLE,  GOLDEN_ANGLE = 2*pi*(1 - 1/PHI)

    The Fibonacci connection here is the appearance of Fibonacci numbers in the
    visible parastichy counts, not in the coordinates themselves.  It is the
    construction usually meant by "Fibonacci distribution" in the phyllotaxis
    literature, and it gives an area-uniform (not logarithmic) radial law.

``"rational_angle"``
    Same radial law as ``"radial"`` but the angular step is the *rational*
    Fibonacci approximant ``theta_step = 2*pi*F_{k-1}/F_k``.  This makes the
    layout exactly ``F_k``-fold periodic in angle, which is the qualitative
    difference between a Fibonacci construction and a golden-ratio one:
    rational steps close, irrational steps never do.

All three are normalised onto the same ``[r_min, r_max]`` so that comparisons
against ``golden_spiral_layout`` hold the radial extent fixed.
"""

from __future__ import annotations

import numpy as np

from .golden_spiral import PHI, _rescale

GOLDEN_ANGLE = 2.0 * np.pi * (1.0 - 1.0 / PHI)  # ~= 2.39996 rad = 137.507 deg

VARIANTS = ("radial", "golden_angle", "rational_angle")


def fibonacci_sequence(n: int, f0: int = 1, f1: int = 1) -> np.ndarray:
    """First ``n`` terms of ``F_k = F_{k-1} + F_{k-2}`` as float64.

    Defaults are ``F_1 = F_2 = 1``.  float64 is exact for Fibonacci numbers up
    to F_78; beyond that the values exceed 2**53 and lose exactness, which does
    not matter after normalisation but is stated here for honesty.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    seq = np.empty(n, dtype=float)
    a, b = float(f0), float(f1)
    for k in range(n):
        seq[k] = a
        a, b = b, a + b
    return seq


def fibonacci_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    variant: str = "radial",
    theta_step: float | None = None,
    theta_0: float = 0.0,
    orientation: float = 0.0,
    skip: int = 0,
    rational_index: int = 8,
) -> np.ndarray:
    """Return an ``(n, 2)`` array of antenna positions in metres.

    Parameters
    ----------
    n : number of antennas.
    r_min, r_max : radial extent in metres (all variants are normalised to it).
    variant : one of :data:`VARIANTS`.
    theta_step : angular increment in radians for ``"radial"``.  ``None`` ->
        :data:`GOLDEN_ANGLE`.
    theta_0 : angle of the first antenna.
    orientation : rigid rotation of the layout.
    skip : number of leading Fibonacci terms to drop before taking ``n`` of them
        (``"radial"``/``"rational_angle"``).  ``skip>0`` removes the repeated
        ``1, 1`` head and the strongly non-asymptotic early terms.
    rational_index : ``k`` in ``theta_step = 2*pi*F_{k-1}/F_k`` for
        ``"rational_angle"``.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if r_max <= r_min:
        raise ValueError("r_max must exceed r_min")
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}, got {variant!r}")

    if variant == "golden_angle":
        idx = np.arange(n, dtype=float)
        r = _rescale(np.sqrt((idx + 0.5) / n), r_min, r_max)
        theta = theta_0 + idx * GOLDEN_ANGLE
    else:
        fib = fibonacci_sequence(n + skip)[skip:]
        r = _rescale(fib, r_min, r_max)
        if variant == "radial":
            step = GOLDEN_ANGLE if theta_step is None else float(theta_step)
        else:
            f = fibonacci_sequence(rational_index + 1)
            step = 2.0 * np.pi * f[rational_index - 1] / f[rational_index]
        theta = theta_0 + np.arange(n, dtype=float) * step

    ang = theta + orientation
    return np.column_stack((r * np.cos(ang), r * np.sin(ang)))

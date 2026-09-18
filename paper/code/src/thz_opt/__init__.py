"""thz_opt -- classical foundation for QUBO-based terahertz interferometric array design.

The package is organised as a pipeline:

    arrays          -> antenna (x, y) layouts in metres
    interferometry  -> baselines, (u, v) coordinates, Earth rotation, PSF
    metrics         -> UV coverage / redundancy / baseline-distribution / PSF metrics
    constraints     -> physical feasibility (separation, PWV placeholder, phase graph)
    qubo            -> binary-selection objective, coefficients, exhaustive validation

Conventions used throughout the package are documented in
``docs/mathematical_formulation.md``.  The short version:

* Antenna coordinates are local East-North-Up (ENU) metres, ``x`` = East, ``y`` = North.
* Baselines are formed for ``i < j`` only, ``b_ij = r_j - r_i`` (N(N-1)/2 of them).
* ``(u, v)`` are in wavelengths.  Conjugate points ``(-u, -v)`` are physically
  sampled as well and are added explicitly where a routine says so.
"""

__version__ = "0.1.0"

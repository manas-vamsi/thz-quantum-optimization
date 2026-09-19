"""A UV weighting *derived* from a science case, not chosen.

The problem this fixes
----------------------
Everywhere else in this project the array is scored by the number of occupied
UV cells, optionally multiplied by a radial weight. That weight was picked by
hand, and optimising a figure of merit you invented yourself is circular: you
can always make an array look good by choosing the weight that likes it.

This module removes the choice. Given a source model and a parameter you want
to measure, estimation theory says exactly how much each UV sample is worth,
and the weighting falls out with nothing left to tune.

The derivation
--------------
For an interferometer the measured visibility at spatial frequency
``(u, v)`` is the Fourier transform of the sky brightness. Suppose the
brightness depends on a parameter ``p`` we want to estimate -- here the depth
of an annular gap in a protoplanetary disc. With independent Gaussian noise of
variance ``sigma^2`` per visibility, the Fisher information about ``p`` from a
set of sampled cells ``S`` is

    F(p) = (1/sigma^2) * sum_{(u,v) in S} |dV(u,v)/dp|^2 ,                (1)

and the Cramer-Rao bound gives ``var(p_hat) >= 1/F(p)``. Maximising the Fisher
information is therefore exactly minimising the best achievable error bar on
``p``. Equation (1) is a *sum over sampled cells of a per-cell weight*, which
is precisely the form the optimiser already maximises -- so the science case
plugs into the existing machinery with no change to the solver, and the weight

    w(u, v) = |dV(u, v)/dp|^2                                             (2)

is derived, not selected. Two arrays can now be compared by the error bar they
deliver on a physical quantity, in the units of that quantity.

The source model
----------------
An axisymmetric disc, so the visibility is the order-zero Hankel transform of
the radial brightness profile,

    V(q) = 2*pi * integral_0^inf I(r) J_0(2*pi*q*r) r dr ,                (3)

with ``r`` in radians on the sky and ``q = |b|/lambda`` in wavelengths. The
profile is a truncated power law with a Gaussian annular depression:

    I(r) = I_0 (r/r_0)^(-gamma) * [1 - p * exp(-(r - r_gap)^2 / (2 s^2))]  (4)

for ``r <= r_out``, zero outside.

Why the derivative must conserve flux
-------------------------------------
Taken literally, ``dI/dp = -I_0 (r/r_0)^(-gamma) * g(r)`` with ``g`` the
Gaussian, and its Hankel transform is largest at ``q = 0``: deepening a gap
removes flux, and nothing measures total flux better than a zero-spacing
sample. The resulting weight peaks at the shortest baseline available, which
is a true statement about a badly posed question -- it says the cheapest way
to detect a gap is to notice the source got fainter.

That is degenerate with the overall brightness scale, which an interferometer
does not measure well and which single-dish photometry supplies independently.
The well-posed question is about *shape at fixed flux*, so the model
renormalises as the gap deepens: ``I(r; p) = A(p) * base(r) * [1 - p g(r)]``
with ``A(p)`` fixed by ``d/dp integral I dOmega = 0``. Differentiating,

    dI/dp = A * base * [ (F_g / F_0) * (1 - p g(r)) - g(r) ]              (5)

where ``F_0 = int base (1 - p g) dOmega`` and ``F_g = int base g dOmega``.
This derivative integrates to zero over the sky, so ``dV/dp`` vanishes at
``q = 0`` by construction and the weight peaks at the spatial frequency of the
ring itself -- which is the scale the array actually has to deliver. The
unconstrained form is kept available as ``flux_conserving=False`` so the two
can be compared rather than asserted.

``p`` enters ``I`` non-linearly once ``A(p)`` is included, but (5) is exact at
the model's own ``p``, which is where the Fisher information is evaluated.

Default parameters
------------------
HL Tau, from the ALMA Partnership long-baseline campaign (ALMA Partnership et
al. 2015, ApJL 808, L3, doi:10.1088/2041-8205/808/1/L3): distance 140 pc, dust
continuum detected to roughly 120 au, prominent dark rings at 13.2, 32.3 and
64.2 au. The 32.3 au gap is the default target. The surface-brightness index
``gamma = 0.5`` and gap width 4.5 au are representative of that paper's fitted
rings; they set the *shape* of the weight, and
:func:`weight_sensitivity_to_profile` exists so a reader can see how much the
conclusion depends on them rather than having to trust the choice.

What is still assumed
---------------------
Noise is taken as independent and identical per visibility, so (1) omits
per-baseline sensitivity differences and atmospheric correlation. The disc is
axisymmetric and face-on. These make the weight a clean function of ``q``
alone; a full treatment would make it two-dimensional and correlated.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import j0

__all__ = [
    "DiskModel",
    "HL_TAU",
    "brightness_profile",
    "gap_derivative_profile",
    "hankel_transform",
    "visibility_profile",
    "fisher_weight_profile",
    "resolution_requirement",
    "uv_cell_weights",
    "fisher_information",
    "depth_error_bar",
    "weight_sensitivity_to_profile",
    "PARAMETERS",
]

#: Gap parameters this module can derive a weighting for.
PARAMETERS = ("depth", "radius", "width")

ARCSEC = np.pi / (180.0 * 3600.0)   # radians per arcsecond


@dataclass(frozen=True)
class DiskModel:
    """An axisymmetric disc with one annular gap.

    Radii are in au, distance in pc; ``au / pc`` is arcseconds by definition,
    which is why no explicit conversion constant appears below.
    """

    distance_pc: float = 140.0
    r_out_au: float = 120.0
    r_in_au: float = 1.0
    gap_radius_au: float = 32.3
    gap_width_au: float = 4.5
    gap_depth: float = 0.5
    gamma: float = 0.5
    r0_au: float = 10.0
    name: str = "HL Tau (ALMA Partnership 2015)"

    def angular(self, r_au: np.ndarray) -> np.ndarray:
        """au -> radians on the sky at this distance."""
        return np.asarray(r_au, dtype=float) / self.distance_pc * ARCSEC

    @property
    def theta_out(self) -> float:
        """Angular radius of the emitting disc, in radians."""
        return float(self.angular(self.r_out_au))

    @property
    def theta_gap(self) -> float:
        """Angular width of the gap feature (FWHM), in radians."""
        return float(self.angular(self.gap_width_au) * 2.3548)


HL_TAU = DiskModel()


# --------------------------------------------------------------------------
# source model, equations (4) and (5)
# --------------------------------------------------------------------------

def _radial_grid(model: DiskModel, n: int = 4096) -> np.ndarray:
    return np.linspace(model.r_in_au, model.r_out_au, n)


def brightness_profile(r_au: np.ndarray, model: DiskModel = HL_TAU) -> np.ndarray:
    """``I(r)`` of equation (4), normalised to unit peak."""
    r = np.asarray(r_au, dtype=float)
    base = (r / model.r0_au) ** (-model.gamma)
    gap = np.exp(-0.5 * ((r - model.gap_radius_au) / model.gap_width_au) ** 2)
    prof = base * (1.0 - model.gap_depth * gap)
    prof = np.where((r >= model.r_in_au) & (r <= model.r_out_au), prof, 0.0)
    peak = prof.max()
    return prof / peak if peak > 0 else prof


def _dh_dtheta(r: np.ndarray, model: DiskModel, parameter: str) -> np.ndarray:
    """``d/dtheta [1 - p g(r)]`` for the parameter being estimated."""
    g = np.exp(-0.5 * ((r - model.gap_radius_au) / model.gap_width_au) ** 2)
    if parameter == "depth":
        return -g
    if parameter == "radius":
        return -model.gap_depth * g * (r - model.gap_radius_au) / model.gap_width_au ** 2
    if parameter == "width":
        return (-model.gap_depth * g * (r - model.gap_radius_au) ** 2
                / model.gap_width_au ** 3)
    raise ValueError(f"unknown parameter {parameter!r}; "
                     f"choose from 'depth', 'radius', 'width'")


def gap_derivative_profile(r_au: np.ndarray, model: DiskModel = HL_TAU,
                           flux_conserving: bool = True,
                           parameter: str = "depth") -> np.ndarray:
    """``dI/dtheta`` of equation (5), for one gap parameter.

    With ``flux_conserving`` the perturbation integrates to zero over the sky,
    so the parameter is measured from the source *shape* and not from its total
    brightness. Set it False to see the degenerate form described above.

    ``parameter`` selects which quantity is being estimated. This matters more
    than it looks: depth, radius and width are informed by different spatial
    frequencies, so "the array for this science case" is only well defined once
    the question is.
    """
    r = np.asarray(r_au, dtype=float)
    inside = (r >= model.r_in_au) & (r <= model.r_out_au)
    base = np.where(inside, (np.maximum(r, model.r_in_au) / model.r0_au)
                    ** (-model.gamma), 0.0)
    dh = _dh_dtheta(r, model, parameter)

    if flux_conserving:
        # solid-angle weighting is r dr for a face-on disc
        rg = _radial_grid(model)
        b = (rg / model.r0_au) ** (-model.gamma)
        g = np.exp(-0.5 * ((rg - model.gap_radius_au) / model.gap_width_au) ** 2)
        h = 1.0 - model.gap_depth * g
        f0 = np.trapezoid(b * h * rg, rg)
        h_theta = np.trapezoid(b * _dh_dtheta(rg, model, parameter) * rg, rg)
        gap_term = 1.0 - model.gap_depth * np.exp(
            -0.5 * ((r - model.gap_radius_au) / model.gap_width_au) ** 2)
        d = base * (dh - (h_theta / f0) * gap_term)
    else:
        d = base * dh

    d = np.where(inside, d, 0.0)
    # divide by the same peak that normalises I(r), so the two are comparable
    rg = _radial_grid(model)
    unnorm = (rg / model.r0_au) ** (-model.gamma) * (
        1.0 - model.gap_depth
        * np.exp(-0.5 * ((rg - model.gap_radius_au) / model.gap_width_au) ** 2))
    scale = unnorm.max() if unnorm.max() > 0 else 1.0
    return d / scale


def hankel_transform(profile: np.ndarray, r_au: np.ndarray, q: np.ndarray,
                     model: DiskModel = HL_TAU) -> np.ndarray:
    """Order-zero Hankel transform, equation (3), for ``q`` in wavelengths.

    Trapezoidal on the supplied radial grid. The integrand is smooth and
    compactly supported, so this converges quickly; the test suite checks it
    against the analytic transform of a uniform disc.
    """
    r_rad = model.angular(r_au)
    qq = np.atleast_1d(np.asarray(q, dtype=float))
    # (n_q, n_r) kernel; fine at the grid sizes used here
    kernel = j0(2.0 * np.pi * qq[:, None] * r_rad[None, :])
    integrand = profile[None, :] * kernel * r_rad[None, :]
    return 2.0 * np.pi * np.trapezoid(integrand, r_rad, axis=1)


def visibility_profile(q: np.ndarray, model: DiskModel = HL_TAU,
                       n_r: int = 4096) -> np.ndarray:
    """``V(q)`` of the disc, equation (3)."""
    r = _radial_grid(model, n_r)
    return hankel_transform(brightness_profile(r, model), r, q, model)


def fisher_weight_profile(q: np.ndarray, model: DiskModel = HL_TAU,
                          n_r: int = 4096, normalise: bool = True,
                          flux_conserving: bool = True,
                          parameter: str = "depth") -> np.ndarray:
    """``w(q) = |dV/dp|^2`` of equation (2).

    This is the whole point of the module: the value of sampling spatial
    frequency ``q``, for the purpose of measuring the gap depth.
    """
    r = _radial_grid(model, n_r)
    dv = hankel_transform(
        gap_derivative_profile(r, model, flux_conserving, parameter), r, q, model)
    w = dv ** 2
    if normalise and w.max() > 0:
        w = w / w.max()
    return w


# --------------------------------------------------------------------------
# what the science demands of the array, independent of any optimiser
# --------------------------------------------------------------------------

def resolution_requirement(model: DiskModel = HL_TAU, wavelength_m: float = 1.0e-3) -> dict:
    """Baseline range the science case implies, before any optimisation.

    The gap must be resolved, which needs a synthesised beam no larger than its
    width, ``lambda / b_max <= theta_gap``; and the disc must not be resolved
    out, which needs a largest recoverable angular scale above the disc
    diameter. Using the standard ``LAS ~ 0.6 lambda / b_min`` convention, both
    give hard baseline limits from the source alone.
    """
    theta_gap = model.theta_gap
    theta_disc = 2.0 * model.theta_out
    b_max = wavelength_m / theta_gap
    b_min = 0.6 * wavelength_m / theta_disc
    return {
        "wavelength_m": wavelength_m,
        "theta_gap_arcsec": theta_gap / ARCSEC,
        "theta_disc_arcsec": theta_disc / ARCSEC,
        "b_max_required_m": b_max,
        "b_min_allowed_m": b_min,
        "q_max_required": b_max / wavelength_m,
        "q_min_allowed": b_min / wavelength_m,
        "source": model.name,
    }


def uv_cell_weights(grid, model: DiskModel = HL_TAU,
                    parameter: str = "depth") -> np.ndarray:
    """Per-cell Fisher weight over a :class:`~thz_opt.interferometry.uv.UVGrid`.

    Returned flattened in the same ``iu * n_cells + iv`` order the occupancy
    arrays use, so it can be handed straight to ``UVState(cell_weights=...)``
    or to the QUBO builders.
    """
    n = grid.n_cells
    edge = (np.arange(n) + 0.5) * (2.0 * grid.uv_max / n) - grid.uv_max
    uu, vv = np.meshgrid(edge, edge, indexing="ij")
    q = np.hypot(uu, vv).ravel()
    return fisher_weight_profile(q, model, parameter=parameter)


def fisher_information(q_sampled: np.ndarray, model: DiskModel = HL_TAU,
                       sigma: float = 1.0, parameter: str = "depth") -> float:
    """``F(theta)`` of equation (1) for a set of sampled spatial frequencies."""
    w = fisher_weight_profile(np.asarray(q_sampled, dtype=float), model,
                              normalise=False, parameter=parameter)
    return float(w.sum() / sigma ** 2)


def depth_error_bar(q_sampled: np.ndarray, model: DiskModel = HL_TAU,
                    sigma: float = 1.0, parameter: str = "depth") -> float:
    """Cramer-Rao lower bound on the gap-depth error, ``1/sqrt(F)``.

    An array is finally comparable in a physical unit: this is the smallest
    error bar on a dimensionless gap depth that any unbiased estimator can
    achieve with the given sampling and noise.
    """
    f = fisher_information(q_sampled, model, sigma, parameter)
    return float("inf") if f <= 0 else float(1.0 / np.sqrt(f))


def weight_sensitivity_to_profile(q: np.ndarray, model: DiskModel = HL_TAU,
                                  gammas=(0.0, 0.5, 1.0),
                                  widths_au=(3.0, 4.5, 6.0)) -> dict:
    """How much the derived weight depends on the two assumed shape parameters.

    The gap radius and the disc size are measured; ``gamma`` and the gap width
    are fitted quantities with real uncertainty. If the weight's peak moved
    substantially across this range the science case would not pin the array
    down, so this is reported rather than hidden.
    """
    out = {}
    for g in gammas:
        for w_au in widths_au:
            m = DiskModel(**{**model.__dict__, "gamma": g, "gap_width_au": w_au})
            w = fisher_weight_profile(q, m)
            out[(g, w_au)] = {
                "q_peak": float(np.asarray(q)[int(np.argmax(w))]),
                "q_half_max_hi": float(
                    np.asarray(q)[int(np.max(np.flatnonzero(w >= 0.5)))]
                    if np.any(w >= 0.5) else np.nan),
            }
    return out

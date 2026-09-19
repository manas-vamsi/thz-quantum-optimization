"""Instrumental effects: primary beam, smearing, and a real noise scale.

What was missing and why it matters
-----------------------------------
The imaging in :mod:`thz_opt.imaging.simulate` treats the interferometer as a
perfect sampler of the sky's Fourier transform with white noise added. Three
omissions from that picture change what an array can actually deliver, and all
three get worse exactly where this project pushes -- long baselines at short
wavelengths.

**Primary beam.** Each antenna sees a patch of sky of angular width about
``lambda / D``, not the whole field. Sources away from the pointing centre are
attenuated, so the usable field of view is set by the dish, not by the
correlator or the grid. For an 8 m dish at 230 GHz the half-power width is 38
arcseconds, smaller than many of the fields a naive simulation would happily
image.

**Bandwidth smearing.** A finite channel width averages visibilities over a
range of ``lambda``, which stretches the response radially outward from the
phase centre in proportion to distance from it. The fractional smear is
``(delta_nu / nu) * theta``, so a 2 % channel blurs a source 20 beams out by
0.4 beams.

**Time smearing.** A finite integration averages over Earth rotation, which
smears tangentially by ``omega_E * tau * theta``. At 10 s integrations that is
7.3e-4 radians of rotation.

Both smearings vanish at the phase centre and grow linearly with distance from
it, which is why they define a field of view rather than a uniform loss, and
why they are implemented here as a position-dependent blur rather than as a
single attenuation factor.

**Noise.** White noise scaled to a chosen image signal-to-noise is fine for
comparing two arrays under identical conditions, but it cannot answer "how long
would this take to observe?". :func:`sefd` and :func:`radiometer_noise` give
the real scale from the measured atmospheric opacity, so integration times can
be quoted in seconds rather than in arbitrary units.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "primary_beam_fwhm_arcsec",
    "primary_beam",
    "apply_primary_beam",
    "smear_image",
    "system_temperature",
    "sefd",
    "radiometer_noise",
    "OMEGA_EARTH",
    "K_BOLTZMANN",
]

ARCSEC = np.pi / (180.0 * 3600.0)
OMEGA_EARTH = 7.2921159e-5      #: rad/s, sidereal rotation rate
K_BOLTZMANN = 1.380649e-23      #: J/K


# --------------------------------------------------------------------------
# primary beam
# --------------------------------------------------------------------------

def primary_beam_fwhm_arcsec(dish_m: float, wavelength_m: float,
                             taper: float = 1.13) -> float:
    """Half-power width of one antenna's response.

    ``taper`` is the usual illumination factor: 1.02 for a uniformly
    illuminated aperture, about 1.13 for the tapered illumination real feeds
    produce. ALMA quotes 1.13 lambda / D, which is the default here.
    """
    return float(taper * wavelength_m / dish_m / ARCSEC)


def primary_beam(grid, dish_m: float, wavelength_m: float,
                 taper: float = 1.13) -> np.ndarray:
    """Gaussian primary-beam attenuation over the image grid.

    A Gaussian is the standard working approximation to the main lobe of an
    Airy pattern; it omits the sidelobes, which matter for wide-field imaging
    and are not modelled here.
    """
    from .sky import image_coordinates

    (ll, mm), _ = image_coordinates(grid)
    fwhm = primary_beam_fwhm_arcsec(dish_m, wavelength_m, taper) * ARCSEC
    sigma = fwhm / 2.3548
    return np.exp(-0.5 * (ll ** 2 + mm ** 2) / sigma ** 2)


def apply_primary_beam(sky, dish_m: float, wavelength_m: float, grid,
                       taper: float = 1.13):
    """Attenuate a sky model by the primary beam, returning a new model."""
    from .sky import SkyModel

    a = primary_beam(grid, dish_m, wavelength_m, taper)
    img = np.asarray(sky.image if hasattr(sky, "image") else sky) * a
    pix = sky.pixel_rad if hasattr(sky, "pixel_rad") else 1.0
    name = getattr(sky, "name", "sky")
    return SkyModel(img, pix, f"{name} x primary beam")


# --------------------------------------------------------------------------
# smearing
# --------------------------------------------------------------------------

def smear_image(image: np.ndarray,
                fractional_bandwidth: float = 0.0,
                integration_s: float = 0.0,
                declination_deg: float = 0.0,
                n_steps: int = 9) -> np.ndarray:
    """Blur an image radially by bandwidth and tangentially by time.

    Both effects are averages over a *geometric transform* of the sky, not
    convolutions with a fixed kernel, and implementing them as such is what
    makes them exact in the no-smearing limit.

    Observing at frequency ``nu`` instead of ``nu_0`` scales every baseline in
    wavelengths by ``nu / nu_0``, which is equivalent to scaling the image by
    ``nu_0 / nu``. Averaging over a rectangular channel is therefore an average
    of the image over scale factors spanning ``1 +/- b/2``. Likewise a finite
    integration averages over the rotation the Earth performs during it. Both
    reduce to a single identity transform when the smear is zero, so this
    routine returns the input exactly rather than leaking a few per cent of the
    flux to interpolation, which an earlier polar-resampling implementation did
    -- it lost up to 24 % of a point source and inflated the total flux by 29 %
    while applying no actual smearing at realistic parameters.

    Scaling changes pixel area by ``s^2``, so each term is divided by it and
    total flux is conserved.
    """
    from scipy.ndimage import map_coordinates

    img = np.asarray(image, dtype=float)
    if fractional_bandwidth <= 0.0 and integration_s <= 0.0:
        return img.copy()

    n = img.shape[0]
    c = (n - 1) / 2.0
    b = max(float(fractional_bandwidth), 0.0)
    dphi = OMEGA_EARTH * np.cos(np.deg2rad(declination_deg)) * max(integration_s, 0.0)

    scales = (np.array([1.0]) if b <= 0
              else np.linspace(1.0 - b / 2.0, 1.0 + b / 2.0, n_steps))
    angles = (np.array([0.0]) if dphi <= 0
              else np.linspace(-dphi / 2.0, dphi / 2.0, n_steps))

    yy, xx = np.mgrid[0:n, 0:n]
    dy, dx = yy - c, xx - c
    out = np.zeros_like(img)
    for s_ in scales:
        for a in angles:
            ca, sa = np.cos(a), np.sin(a)
            # where each output pixel came from in the unsmeared sky
            sy = (dy * ca + dx * sa) / s_ + c
            sx = (-dy * sa + dx * ca) / s_ + c
            out += map_coordinates(img, [sy.ravel(), sx.ravel()], order=1,
                                   mode="constant", cval=0.0
                                   ).reshape(img.shape) / s_ ** 2
    return out / (len(scales) * len(angles))


# --------------------------------------------------------------------------
# noise, on a real scale
# --------------------------------------------------------------------------

def system_temperature(tau_zenith: float, elevation_deg: float = 60.0,
                       t_atmosphere_k: float = 265.0,
                       t_receiver_k: float = 60.0,
                       t_cmb_k: float = 2.73) -> dict:
    """System temperature above the atmosphere, from the zenith opacity.

    ``tau_zenith`` is the measured zenith optical depth -- at Hanle a 220 GHz
    radiometer has recorded it since 1999, so this is a measured input rather
    than a guess. The airmass is ``1 / sin(elevation)``, the plane-parallel
    approximation, which is adequate above about 20 degrees.
    """
    if not 0.0 < elevation_deg <= 90.0:
        raise ValueError("elevation must lie in (0, 90] degrees")
    airmass = 1.0 / np.sin(np.deg2rad(elevation_deg))
    tau = float(tau_zenith) * airmass
    transmission = float(np.exp(-tau))
    t_sky = t_atmosphere_k * (1.0 - transmission) + t_cmb_k * transmission
    t_sys_above = (t_receiver_k + t_sky) / transmission
    return {
        "tau_zenith": float(tau_zenith),
        "airmass": float(airmass),
        "tau_line_of_sight": tau,
        "transmission": transmission,
        "t_sky_k": float(t_sky),
        "t_sys_k": float(t_sys_above),
    }


def sefd(t_sys_k: float, dish_m: float, aperture_efficiency: float = 0.7) -> float:
    """System equivalent flux density in janskys.

    ``SEFD = 2 k T_sys / A_eff`` with ``A_eff = eta * pi D^2 / 4``.
    """
    area = aperture_efficiency * np.pi * dish_m ** 2 / 4.0
    return float(2.0 * K_BOLTZMANN * t_sys_k / area / 1e-26)


def radiometer_noise(sefd_jy: float, n_antennas: int, bandwidth_hz: float,
                     integration_s: float, n_polarizations: int = 2) -> dict:
    """Point-source sensitivity of a whole array, in janskys.

    ``sigma = SEFD / (eta_c sqrt(n_pol N (N-1) delta_nu tau))``. The
    correlator efficiency is taken as 1, which is optimistic by a few per cent.
    """
    if n_antennas < 2:
        raise ValueError("an interferometer needs at least two antennas")
    n_baselines = n_antennas * (n_antennas - 1) / 2
    denom = np.sqrt(n_polarizations * n_baselines * bandwidth_hz * integration_s)
    sigma = float(sefd_jy / denom) if denom > 0 else float("inf")
    return {
        "sefd_jy": float(sefd_jy),
        "n_baselines": int(n_baselines),
        "sigma_jy": sigma,
        "integration_s": float(integration_s),
        "bandwidth_ghz": bandwidth_hz / 1e9,
    }

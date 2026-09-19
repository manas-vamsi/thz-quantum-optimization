"""End-to-end imaging: make a picture, then measure how wrong it is.

Everything else in this project scores an array by a *predictor* of image
quality -- occupied UV cells, sidelobe energy, a Fisher bound. None of them is
an image. The whole framework rests on the assumption that improving those
numbers improves the picture, and that assumption has never been tested here.

This package closes the loop:

    sky model -> visibilities -> sample on the array's UV coverage
              -> add noise -> dirty image -> CLEAN -> compare to the truth

so two arrays can be ranked by the fidelity of what they actually reconstruct.
If the ranking agrees with the UV-cell ranking, the cheap metric is validated.
If it disagrees, the cheap metric was misleading and that is worth more than
another coverage table.

The transform convention matches
:func:`thz_opt.interferometry.psf.dirty_beam` exactly, so the dirty beam used
for deconvolution here is the same object reported elsewhere.
"""

from .sky import (
    SkyModel,
    disc_with_gap_image,
    gaussian_source,
    point_sources,
    pixel_scale_arcsec,
)
from .simulate import observe, predict_visibilities, thermal_noise_sigma
from .clean import hogbom_clean, restore
from .fidelity import fidelity_metrics
from .instrument import (
    apply_primary_beam,
    primary_beam,
    primary_beam_fwhm_arcsec,
    radiometer_noise,
    sefd,
    smear_image,
    system_temperature,
)

__all__ = [
    "SkyModel",
    "point_sources",
    "gaussian_source",
    "disc_with_gap_image",
    "pixel_scale_arcsec",
    "predict_visibilities",
    "thermal_noise_sigma",
    "observe",
    "hogbom_clean",
    "restore",
    "fidelity_metrics",
    "primary_beam",
    "primary_beam_fwhm_arcsec",
    "apply_primary_beam",
    "smear_image",
    "system_temperature",
    "sefd",
    "radiometer_noise",
]

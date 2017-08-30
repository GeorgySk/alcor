import enum
from typing import List

import numpy as np

from alcor.models import Star
from alcor.models.star import SpectralTypeEnum


def assign_magnitudes(stars: List[Star],
                      max_carbon_oxygen_core_wd_mass: float = 1.14,
                      db_to_da_fraction: float = 0.2) -> None:
    co_wds = [star
              for star in stars
              if star.mass < max_carbon_oxygen_core_wd_mass]
    one_wds = [star
               for star in stars
               if star.mass >= max_carbon_oxygen_core_wd_mass]

    for star in co_wds:
        if get_spectral_type(db_to_da_fraction) == SpectralTypeEnum.DA:
            star.spectral_type = SpectralTypeEnum.DA
            # da_interpolation()
        else:
            star.spectral_type = SpectralTypeEnum.DB
            # db_interpolation()

    for star in one_wds:
        star.spectral_type = SpectralTypeEnum.ONe
        # one_interpolation()


def get_spectral_type(db_to_da_fraction: float) -> enum.Enum:
    if np.random.rand() < db_to_da_fraction:
        return SpectralTypeEnum.DB
    return SpectralTypeEnum.DA

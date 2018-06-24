from functools import partial

import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

from alcor.models.star import GalacticDiskType
from alcor.types import ArrayOperatorType
from .utils import (immutable_array,
                    linear_function)

# TODO: we have this in another module. take to utils after merge
linear_estimation = partial(InterpolatedUnivariateSpline,
                            k=1)

# values from Althaus priv. comm (X = 0.725, Y = 0.265)
MODEL_SOLAR_MASSES = immutable_array([1.00, 1.50, 1.75, 2.00, 2.25,
                                      2.50, 3.00, 3.50, 4.00, 5.00])
MODEL_SOLAR_TIMES = immutable_array([8.614, 1.968, 1.249, 0.865, 0.632,
                                     0.480, 0.302, 0.226, 0.149, 0.088])
# values from Althaus priv. comm (X = 0.752, Y = 0.247)
MODEL_SUBSOLAR_MASSES = immutable_array([0.85, 1.00, 1.25, 1.50,
                                         1.75, 2.00, 3.00])
MODEL_SUBSOLAR_TIMES = immutable_array([10.34, 5.756, 2.623, 1.412,
                                        0.905, 0.639, 0.245])
LOW_MASS_FUNCTION = partial(linear_function,
                            factor=0.096,
                            const=0.429)
MEDIUM_MASS_FUNCTION = partial(linear_function,
                               factor=0.137,
                               const=0.3183)
HIGH_MASS_FUNCTION = partial(linear_function,
                             factor=0.1057,
                             const=0.5061)


def get_white_dwarfs(stars: pd.DataFrame,
                     *,
                     max_galactic_structure_age: float,
                     chandrasekhar_limit: float,
                     solar_metallicity: float,
                     subsolar_metallicity: float,
                     min_cooling_time: float = 0.) -> pd.DataFrame:
    """
    Filters white dwarfs stars from initial sample of main sequence stars
    and assigns metallicities, cooling times and masses.

    :param stars: main sequence stars
    :param max_galactic_structure_age: the highest age of thin disk,
    thick disk and halo
    :param chandrasekhar_limit: maximum mass of a stable white dwarf
    :param solar_metallicity: metallicity assigned to all thin
    and thick disks white dwarfs due to relatively young ages
    :param subsolar_metallicity: metallicity assigned to all halo white dwarfs
    :param min_cooling_time: natural lower limit for cooling time
    :return: white dwarfs
    """
    stars['metallicity'] = get_metallicities(
            galactic_disks_types=stars['galactic_disk_type'].values,
            subsolar_metallicity=subsolar_metallicity,
            solar_metallicity=solar_metallicity)

    halo_stars_mask = stars['galactic_disk_type'] == GalacticDiskType.halo
    disks_stars_mask = ~halo_stars_mask

    stars.loc[disks_stars_mask, 'lifetime'] = main_sequence_lifetimes(
            stars.loc[disks_stars_mask, 'progenitor_mass'].values,
            model_masses=MODEL_SOLAR_MASSES,
            model_times=MODEL_SOLAR_TIMES)
    stars.loc[halo_stars_mask, 'lifetime'] = main_sequence_lifetimes(
            stars.loc[halo_stars_mask, 'progenitor_mass'].values,
            model_masses=MODEL_SUBSOLAR_MASSES,
            model_times=MODEL_SUBSOLAR_TIMES)

    stars['cooling_time'] = (max_galactic_structure_age
                             - stars['birth_time']
                             - stars['lifetime'])
    stars.drop(columns='lifetime',
               inplace=True)

    stars = stars[stars['cooling_time'] > min_cooling_time]

    stars['mass'] = white_dwarf_masses(stars['progenitor_mass'].values)

    return stars[stars['mass'] <= chandrasekhar_limit]


def get_metallicities(*,
                      galactic_disks_types: np.ndarray,
                      subsolar_metallicity: float,
                      solar_metallicity: float) -> np.ndarray:
    result = np.empty(galactic_disks_types.size)

    halo_mask = np.equal(galactic_disks_types.astype(int),
                         np.full(galactic_disks_types.size,
                                 GalacticDiskType.halo.value))
    # halo_mask = galactic_disks_types == GalacticDiskType.halo.value

    result[halo_mask] = subsolar_metallicity
    result[~halo_mask] = solar_metallicity

    return result


def main_sequence_lifetimes(masses: np.ndarray,
                            *,
                            model_masses: np.ndarray,
                            model_times: np.ndarray) -> np.ndarray:
    """
    Calculates lifetime of a main sequence star
    according to model by Leandro & Renedo et al.(2010).
    """
    if masses.size == 0:
        return np.array([])

    spline = linear_estimation(x=model_masses,
                               y=model_times)

    return estimated_times(
            masses=masses,
            spline=spline,
            rightmost_model_mass=model_masses[-1],
            rightmost_model_time=model_times[-1])


# TODO: this looks like 'white_dwarf_masses' function
def estimated_times(*,
                    masses: np.ndarray,
                    spline: InterpolatedUnivariateSpline,
                    rightmost_model_mass: float,
                    rightmost_model_time: float) -> np.ndarray:
    extrapolating_function = partial(extrapolate_main_sequence_lifetimes,
                                     rightmost_mass=rightmost_model_mass,
                                     rightmost_time=rightmost_model_time)

    lt_max_masses_mask = masses < rightmost_model_mass

    return np.piecewise(masses,
                        condlist=[lt_max_masses_mask, ~lt_max_masses_mask],
                        funclist=[spline, extrapolating_function])


def extrapolate_main_sequence_lifetimes(masses: np.ndarray,
                                        *,
                                        rightmost_mass: float,
                                        rightmost_time: float) -> np.ndarray:
    """
    Extrapolates main sequence stars' (progenitors)
    lifetime vs mass to the right.
    Unlike linear extrapolation this one makes sure
    that no negative values will be produced
    by the fact that 1 / x > 0 for x → ∞
    """
    return rightmost_time * rightmost_mass / masses


def white_dwarf_masses(
        progenitor_masses: np.ndarray,
        *,
        low_mass: float = 2.7,
        high_mass: float = 6.,
        low_mass_function: ArrayOperatorType = LOW_MASS_FUNCTION,
        medium_mass_function: ArrayOperatorType = MEDIUM_MASS_FUNCTION,
        high_mass_function: ArrayOperatorType = HIGH_MASS_FUNCTION
        ) -> np.ndarray:
    """
    IFMR (Initial-to-Final Mass Relation)
    according to model by Catalan et al. 2008.
    Combination with a model by Iben for masses greater than 6 solar masses
    """
    low_progenitor_masses_mask = progenitor_masses < low_mass
    medium_progenitor_masses_mask = ((progenitor_masses >= low_mass)
                                     & (progenitor_masses <= high_mass))
    high_progenitor_masses_mask = (progenitor_masses > high_mass)

    masks = [low_progenitor_masses_mask,
             medium_progenitor_masses_mask,
             high_progenitor_masses_mask]
    functions = [low_mass_function,
                 medium_mass_function,
                 high_mass_function]

    return np.piecewise(progenitor_masses,
                        condlist=masks,
                        funclist=functions)

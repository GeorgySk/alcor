from functools import partial
from math import log10
from typing import (Dict,
                    Tuple,
                    List)

import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

from alcor.models.star import SpectralType

GRAVITATIONAL_CONST_CM_S_KG = 6.67e-5
SOLAR_MASS_KG = 1.989e30

linear_estimation = partial(InterpolatedUnivariateSpline,
                            k=1)


def assign_estimated_values(
        stars: pd.DataFrame,
        *,
        max_carbon_oxygen_core_wd_mass: float = 1.14,
        db_to_da_fraction: float = 0.2,
        da_cooling_sequences: Dict[int, Dict[str, np.ndarray]],
        da_color_table: Dict[str, np.ndarray],
        db_cooling_sequences: Dict[int, Dict[str, np.ndarray]],
        db_color_table: Dict[str, np.ndarray],
        one_color_table: Dict[str, np.ndarray]
        # TODO: we should return DataFrame
        ) -> List[pd.Series]:
    cooling_sequences = {SpectralType.DA: da_cooling_sequences,
                         SpectralType.DB: db_cooling_sequences}
    colo_tables = {SpectralType.DA: da_color_table,
                   SpectralType.DB: db_color_table}
    metallicities = {SpectralType.DA: [0.001, 0.01, 0.03, 0.06],
                     SpectralType.DB: [0.001, 0.01, 0.06]}

    carbon_oxygen_white_dwarfs_mask = (stars['mass']
                                       < max_carbon_oxygen_core_wd_mass)
    carbon_oxygen_white_dwarfs = stars[carbon_oxygen_white_dwarfs_mask]
    oxygen_neon_white_dwarfs = stars[~carbon_oxygen_white_dwarfs_mask]

    for _, star in carbon_oxygen_white_dwarfs.iterrows():
        spectral_type = generate_spectral_type(db_to_da_fraction)
        star['spectral_type'] = spectral_type

        set_estimations_to_da_db_white_dwarf(
                star,
                cooling_sequences=cooling_sequences[spectral_type],
                color_table=colo_tables[spectral_type],
                metallicities=metallicities[spectral_type])

    oxygen_neon_white_dwarfs['spectral_type'] = SpectralType.ONe

    for _, star in oxygen_neon_white_dwarfs.iterrows():
        set_estimations_to_oxygen_neon_white_dwarf(star,
                                                   color_table=one_color_table)

    return carbon_oxygen_white_dwarfs + oxygen_neon_white_dwarfs


def generate_spectral_type(db_to_da_fraction: float) -> SpectralType:
    if np.random.rand() < db_to_da_fraction:
        return SpectralType.DB
    return SpectralType.DA


def set_estimations_to_oxygen_neon_white_dwarf(
        star: pd.Series,
        *,
        color_table: Dict[str, np.ndarray],
        one_model: bool = True,
        by_logarithm: bool = False) -> None:
    star_mass = star['mass']
    star_cooling_time = star['cooling_time']
    mass_grid = color_table['mass']
    cooling_time_grid = color_table['cooling_time']
    pre_wd_lifetime_grid = color_table['pre_wd_lifetime_grid']

    min_mass_index = calculate_index(star_mass,
                                     grid=mass_grid)

    if star_mass < mass_grid[0] or star_mass >= mass_grid[-1]:
        estimate_interest_value = extrapolate_interest_value
    else:
        estimate_interest_value = interpolate_interest_value

    estimate = partial(estimate_interest_value,
                       star_mass=star_mass,
                       star_cooling_time=star_cooling_time,
                       min_mass=mass_grid[min_mass_index],
                       max_mass=mass_grid[min_mass_index + 1],
                       min_mass_index=min_mass_index,
                       cooling_time_grid=cooling_time_grid,
                       pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                       by_logarithm=by_logarithm,
                       one_model=one_model)

    star['luminosity'] = estimate(interest_sequence_grid='luminosity')
    v_ubvri_absolute = estimate(interest_sequence_grid='v_ubvri_absolute')
    bv_ubvri = estimate(interest_sequence_grid='bv_ubvri')
    vi_ubvri = estimate(interest_sequence_grid='vi_ubvri')
    vr_ubvri = estimate(interest_sequence_grid='vr_ubvri')
    uv_ubvri = estimate(interest_sequence_grid='uv_ubvri')
    log_effective_temperature = estimate(
            interest_sequence_grid='log_effective_temperature')

    star['effective_temperature'] = 10. ** log_effective_temperature

    star['u_ubvri_absolute'] = uv_ubvri + v_ubvri_absolute
    star['b_ubvri_absolute'] = bv_ubvri + v_ubvri_absolute
    star['r_ubvri_absolute'] = v_ubvri_absolute - vr_ubvri
    star['i_ubvri_absolute'] = v_ubvri_absolute - vi_ubvri
    star['v_ubvri_absolute'] = v_ubvri_absolute


def set_estimations_to_da_db_white_dwarf(
        star: pd.Series,
        *,
        cooling_sequences: Dict[int, Dict[str, np.ndarray]],
        color_table: Dict[str, np.ndarray],
        metallicities: List[float]) -> pd.Series:
    star_metallicity = star['metallicity']
    star_mass = star['mass']
    star_cooling_time = star['cooling_time']

    min_metallicity_index = get_min_metallicity_index(
            star_metallicity=star_metallicity,
            grid_metallicities=metallicities)
    min_metallicity = metallicities[min_metallicity_index]
    max_metallicity = metallicities[min_metallicity_index + 1]

    min_metallicity_grids = cooling_sequences[int(
            min_metallicity * 1e3)]
    max_metallicity_grids = cooling_sequences[int(
            max_metallicity * 1e3)]

    estimate = partial(estimate_edge_case,
                       star_mass=star_mass,
                       star_cooling_time=star_cooling_time)

    min_luminosity = estimate(
            cooling_sequences=min_metallicity_grids,
            interest_sequence_str='luminosity',
            by_logarithm=False)
    max_luminosity = estimate(
            cooling_sequences=max_metallicity_grids,
            interest_sequence_str='luminosity',
            by_logarithm=False)
    min_effective_temperature = estimate(
            cooling_sequences=min_metallicity_grids,
            interest_sequence_str='effective_temperature',
            by_logarithm=True)
    max_effective_temperature = estimate(
            cooling_sequences=max_metallicity_grids,
            interest_sequence_str='effective_temperature',
            by_logarithm=True)

    star['luminosity'] = -estimate_at(star_metallicity,
                                      x=(min_metallicity, max_metallicity),
                                      y=(min_luminosity, max_luminosity))

    star['effective_temperature'] = estimate_at(
            star_metallicity,
            x=(min_metallicity, max_metallicity),
            y=(min_effective_temperature, max_effective_temperature))

    return star_with_colors(star,
                            color_table=color_table)


def estimate_edge_case(*,
                       cooling_sequences: Dict[str, np.ndarray],
                       star_mass: float,
                       star_cooling_time: float,
                       interest_sequence_str: str,
                       by_logarithm: bool) -> float:
    mass_grid = cooling_sequences['mass']
    cooling_time_grid = cooling_sequences['cooling_time']
    pre_wd_lifetime_grid = cooling_sequences['pre_wd_lifetime_grid']
    mass_index = calculate_index(star_mass,
                                 grid=mass_grid)

    if star_mass < mass_grid[0] or star_mass >= mass_grid[-1]:
        estimate_interest_value = extrapolate_interest_value
    else:
        estimate_interest_value = interpolate_interest_value

    return estimate_interest_value(
            star_mass=star_mass,
            star_cooling_time=star_cooling_time,
            cooling_time_grid=cooling_time_grid,
            pre_wd_lifetime_grid=pre_wd_lifetime_grid,
            min_mass_index=mass_index,
            interest_sequence_grid=cooling_sequences[interest_sequence_str],
            min_mass=mass_grid[mass_index],
            max_mass=mass_grid[mass_index + 1],
            by_logarithm=by_logarithm)


def star_with_colors(star: pd.Series,
                     *,
                     color_table: Dict[str, np.ndarray]) -> pd.Series:
    star_mass = star['mass']
    star_luminosity = star['luminosity']
    luminosity_grid = color_table['luminosity']
    mass_grid = color_table['mass_grid']

    min_mass_index = calculate_index(star_mass,
                                     grid=mass_grid)
    max_mass_index = min_mass_index + 1

    colors = ['u_ubvri_absolute',
              'b_ubvri_absolute',
              'v_ubvri_absolute',
              'r_ubvri_absolute',
              'i_ubvri_absolute']

    min_luminosity_grid = luminosity_grid[min_mass_index, :]
    max_luminosity_grid = luminosity_grid[max_mass_index, :]

    row_index = calculate_index(star_luminosity,
                                grid=min_luminosity_grid)
    next_row_index = calculate_index(
            star_luminosity,
            grid=max_luminosity_grid)

    if (star_luminosity > min_luminosity_grid[0] or
            star_luminosity > max_luminosity_grid[0]):
        min_mass = mass_grid[min_mass_index]
        max_mass = mass_grid[max_mass_index]
    elif (star_luminosity < min_luminosity_grid[-1] or
            star_luminosity < max_luminosity_grid[-1]):
        min_mass = mass_grid[min_mass_index]
        max_mass = mass_grid[max_mass_index]
    else:
        # TODO: check these indexes, they look suspicious
        min_mass = mass_grid[0]
        max_mass = mass_grid[1]

    for color in colors:
        magnitude_grid = color_table[color]
        min_magnitude_grid = magnitude_grid[min_mass_index, :]
        max_magnitude_grid = magnitude_grid[max_mass_index, :]

        min_magnitude = estimate_at(
                star_luminosity,
                x=(min_luminosity_grid[row_index],
                   min_luminosity_grid[row_index + 1]),
                y=(min_magnitude_grid[row_index],
                   magnitude_grid[row_index + 1]))
        max_magnitude = estimate_at(
                star_luminosity,
                x=(max_luminosity_grid[next_row_index],
                   max_luminosity_grid[next_row_index + 1]),
                y=(max_magnitude_grid[next_row_index],
                   max_magnitude_grid[next_row_index + 1]))
        magnitude = estimate_at(star_mass,
                                x=(min_mass, max_mass),
                                y=(min_magnitude, max_magnitude))
        star[color] = max(0., magnitude)

    return star


def extrapolate_interest_value(*,
                               star_mass: float,
                               star_cooling_time: float,
                               cooling_time_grid: np.ndarray,
                               pre_wd_lifetime_grid: np.ndarray,
                               interest_sequence_grid: np.ndarray,
                               min_mass_index: int,
                               min_mass: float,
                               max_mass: float,
                               by_logarithm: bool,
                               one_model: bool = False) -> float:
    interest_value = partial(get_interest_value,
                             star_cooling_time=star_cooling_time,
                             by_logarithm=by_logarithm,
                             one_model=one_model)

    max_mass_index = min_mass_index + 1

    min_interest_value = interest_value(
            mass_index=min_mass_index,
            cooling_time_grid=cooling_time_grid[min_mass_index, :],
            interest_sequence_grid=interest_sequence_grid[min_mass_index, :],
            pre_wd_lifetime=pre_wd_lifetime_grid[min_mass_index])
    max_interest_value = interest_value(
            mass_index=max_mass_index,
            cooling_time_grid=cooling_time_grid[max_mass_index, :],
            interest_sequence_grid=interest_sequence_grid[max_mass_index, :],
            pre_wd_lifetime=pre_wd_lifetime_grid[max_mass_index])
    return estimate_at(star_mass,
                       x=(min_mass, max_mass),
                       y=(min_interest_value, max_interest_value))


def interpolate_interest_value(*,
                               star_mass: float,
                               star_cooling_time: float,
                               min_mass: float,
                               max_mass: float,
                               min_mass_index: int,
                               cooling_time_grid: np.ndarray,
                               pre_wd_lifetime_grid: np.ndarray,
                               interest_sequence_grid: np.ndarray,
                               by_logarithm: bool,
                               one_model: bool = False) -> float:
    max_mass_index = min_mass_index + 1

    min_cooling_time_grid = cooling_time_grid[min_mass_index, :]
    max_cooling_time_grid = cooling_time_grid[max_mass_index, :]
    min_pre_wd_lifetime = pre_wd_lifetime_grid[min_mass_index]
    max_pre_wd_lifetime = pre_wd_lifetime_grid[max_mass_index]
    min_interest_sequence = interest_sequence_grid[min_mass_index, :]
    max_interest_sequence = interest_sequence_grid[max_mass_index, :]

    if one_model:
        min_extrapolated_interest_value = partial(
                estimated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=(min_cooling_time_grid
                                   + min_pre_wd_lifetime),
                interest_sequence_grid=min_interest_sequence)
        max_extrapolated_interest_value = partial(
                estimated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=(max_cooling_time_grid
                                   + max_pre_wd_lifetime),
                interest_sequence_grid=max_interest_sequence)
    elif by_logarithm:
        min_extrapolated_interest_value = partial(
                estimated_log_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=log10(min_cooling_time_grid
                                        + min_pre_wd_lifetime),
                interest_sequence_grid=log10(min_interest_sequence))
        max_extrapolated_interest_value = partial(
                estimated_log_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=log10(max_cooling_time_grid
                                        + max_pre_wd_lifetime),
                interest_sequence_grid=log10(max_interest_sequence))
    else:
        min_extrapolated_interest_value = partial(
                estimated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=log10(min_cooling_time_grid
                                        + min_pre_wd_lifetime),
                interest_sequence_grid=min_interest_sequence)
        max_extrapolated_interest_value = partial(
                estimated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=log10(max_cooling_time_grid
                                        + max_pre_wd_lifetime),
                interest_sequence_grid=max_interest_sequence)

    extrapolating_by_min_cooling_time_grid = extrapolating_by_grid(
            star_cooling_time,
            cooling_time_grid=min_cooling_time_grid)
    min_row_index = calculate_index(star_cooling_time,
                                    grid=min_cooling_time_grid)

    if extrapolating_by_min_cooling_time_grid:
        x_1 = min_extrapolated_interest_value(min_row_index=min_row_index)
    else:
        y_1 = min_cooling_time_grid[min_row_index]
        y_2 = min_cooling_time_grid[min_row_index + 1]
        x_1 = min_interest_sequence[min_row_index]
        x_2 = min_interest_sequence[min_row_index + 1]

    extrapolating_by_max_cooling_time_grid = extrapolating_by_grid(
            star_cooling_time,
            cooling_time_grid=max_cooling_time_grid)
    max_row_index = calculate_index(star_cooling_time,
                                    grid=max_cooling_time_grid)
    if extrapolating_by_max_cooling_time_grid:
        x_3 = max_extrapolated_interest_value(min_row_index=max_row_index)
    else:
        y_3 = min_cooling_time_grid[max_row_index]
        y_4 = min_cooling_time_grid[max_row_index + 1]
        x_3 = min_interest_sequence[max_row_index]
        x_4 = min_interest_sequence[max_row_index + 1]

    if (not extrapolating_by_min_cooling_time_grid and
            not extrapolating_by_max_cooling_time_grid):
        ym_1 = estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(y_1, y_3))
        ym_2 = estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(y_2, y_4))
        xm_1 = estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(x_1, x_3))
        xm_2 = estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(x_2, x_4))

        return estimate_at(star_cooling_time,
                           x=(ym_1, ym_2),
                           y=(xm_1, xm_2))

    if (not extrapolating_by_min_cooling_time_grid
            and extrapolating_by_max_cooling_time_grid):
        xm_1 = estimate_at(star_cooling_time,
                           x=(y_1, y_2),
                           y=(x_1, x_2))
        return estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(xm_1, x_3))

    if (extrapolating_by_min_cooling_time_grid and
            not extrapolating_by_max_cooling_time_grid):
        xm_2 = estimate_at(star_cooling_time,
                           x=(y_3, y_4),
                           y=(x_3, x_4))
        return estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(x_1, xm_2))

    return estimate_at(star_mass,
                       x=(min_mass, max_mass),
                       y=(x_1, x_3))


def extrapolating_by_grid(star_cooling_time: float,
                          cooling_time_grid: np.ndarray) -> bool:
    if (star_cooling_time < cooling_time_grid[0] or
            star_cooling_time >= cooling_time_grid[-1]):
        return True
    else:
        return False


def get_interest_value(*,
                       star_cooling_time: float,
                       cooling_time_grid: np.ndarray,
                       pre_wd_lifetime: float,
                       interest_sequence_grid: np.ndarray,
                       by_logarithm: bool,
                       one_model: bool = False) -> float:
    min_row_index = calculate_index(star_cooling_time,
                                    grid=cooling_time_grid)

    if (star_cooling_time < cooling_time_grid[0] or
            star_cooling_time > cooling_time_grid[-1]):
        if one_model:
            cooling_time_grid = cooling_time_grid + pre_wd_lifetime
        elif by_logarithm:
            cooling_time_grid = np.log10(cooling_time_grid + pre_wd_lifetime)
            interest_sequence_grid = np.log10(interest_sequence_grid)
        else:
            cooling_time_grid = np.log10(cooling_time_grid + pre_wd_lifetime)

        estimated_value = estimated_interest_value(
                min_row_index=min_row_index,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                interest_sequence_grid=interest_sequence_grid)

        if by_logarithm:
            return 10. ** estimated_value
        else:
            return estimated_value
    else:
        return estimated_interest_value(
                min_row_index=min_row_index,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                interest_sequence_grid=interest_sequence_grid)


def estimated_interest_value(*,
                             star_cooling_time: float,
                             cooling_time_grid: np.ndarray,
                             interest_sequence_grid: np.ndarray,
                             min_row_index: int) -> float:
    return estimate_at(star_cooling_time,
                       x=(cooling_time_grid[min_row_index],
                          cooling_time_grid[min_row_index + 1]),
                       y=(interest_sequence_grid[min_row_index],
                          interest_sequence_grid[min_row_index + 1]))


def estimated_log_interest_value(
        *,
        star_cooling_time: float,
        cooling_time_grid: np.ndarray,
        interest_sequence_grid: np.ndarray,
        min_row_index: int) -> float:
    return 10. ** estimated_interest_value(
            star_cooling_time=star_cooling_time,
            cooling_time_grid=cooling_time_grid,
            interest_sequence_grid=interest_sequence_grid,
            min_row_index=min_row_index)


def get_min_metallicity_index(*,
                              star_metallicity: float,
                              grid_metallicities: List[float]) -> int:
    if (star_metallicity < grid_metallicities[0] or
            star_metallicity > grid_metallicities[-1]):
        raise ValueError(f'There is no support for metallicities '
                         f'lying out of the range of {grid_metallicities}')
    star_metallicity = np.array([star_metallicity])
    left_index = np.searchsorted(grid_metallicities, star_metallicity) - 1.
    return np.asscalar(left_index)


def calculate_index(value: float,
                    *,
                    grid: np.ndarray) -> int:
    if value <= grid[0]:
        return 0
    elif value > grid[-1]:
        # Index of element before the last one
        return -2
    star_cooling_time = np.array([value])
    left_index = np.searchsorted(grid, star_cooling_time) - 1
    return np.asscalar(left_index)


def estimate_at(x_0: float,
                *,
                x: Tuple[float, float],
                y: Tuple[float, float]) -> float:
    spline = linear_estimation(x=x,
                               y=y)
    return spline(x_0)

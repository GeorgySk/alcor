from functools import partial
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
        da_cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        da_color_table: Dict[int, pd.DataFrame],
        db_cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        db_color_table: Dict[int, pd.DataFrame],
        one_color_table: Dict[int, pd.DataFrame]
        # TODO: we should return DataFrame
        ) -> List[pd.Series]:
    # TODO: probably I should save it as metadata inside hdf5 files
    da_pre_wd_lifetimes = {1: np.zeros(7),
                           10: np.zeros(10),
                           30: np.array([11.117, 2.7004, 1.699, 1.2114,
                                         0.9892, 0.7422, 0.4431, 0.0]),
                           60: np.array([11.117, 2.7004, 1.699, 1.2114,
                                         0.9892, 0.7422, 0.4431, 0.0])}
    db_pre_wd_lifetimes = {1: np.zeros(7),
                           10: np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892,
                                         0.7422, 0.4431, 0.287, 0.114]),
                           60: np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892,
                                         0.7422, 0.4431, 0.287, 0.114])}
    pre_wd_lifetimes = {SpectralType.DA: da_pre_wd_lifetimes,
                        SpectralType: db_pre_wd_lifetimes}
    cooling_sequences = {SpectralType.DA: da_cooling_sequences,
                         SpectralType.DB: db_cooling_sequences}
    color_tables = {SpectralType.DA: da_color_table,
                    SpectralType.DB: db_color_table}
    da_int_metallicities = sorted(list(da_cooling_sequences.keys()))
    db_int_metallicities = sorted(list(db_cooling_sequences.keys()))
    da_metallicities = [metallicity / 1e3
                        for metallicity in da_int_metallicities]
    db_metallicities = [metallicity / 1e3
                        for metallicity in db_int_metallicities]
    metallicities = {SpectralType.DA: da_metallicities,
                     SpectralType.DB: db_metallicities}

    carbon_oxygen_white_dwarfs_mask = (stars['mass']
                                       < max_carbon_oxygen_core_wd_mass)
    carbon_oxygen_white_dwarfs = stars[carbon_oxygen_white_dwarfs_mask]
    oxygen_neon_white_dwarfs = stars[~carbon_oxygen_white_dwarfs_mask]

    carbon_oxygen_white_dwarfs['spectral_type'] = generate_spectral_types(
            db_to_da_fraction=db_to_da_fraction,
            size=carbon_oxygen_white_dwarfs.shape[0])

    da_white_dwarfs_mask = (carbon_oxygen_white_dwarfs['spectral_type']
                            == SpectralType.DA)
    da_white_dwarfs = carbon_oxygen_white_dwarfs[da_white_dwarfs_mask]
    db_white_dwarfs = carbon_oxygen_white_dwarfs[~da_white_dwarfs_mask]

    white_dwarfs_by_spectral_types = {SpectralType.DA: da_white_dwarfs,
                                      SpectralType.DB: db_white_dwarfs}

    colors = ['u_ubvri_absolute',
              'b_ubvri_absolute',
              'v_ubvri_absolute',
              'r_ubvri_absolute',
              'i_ubvri_absolute',
              'j_ubvri_absolute']

    for spectral_type, white_dwarfs in white_dwarfs_by_spectral_types.items():
        # TODO: why is there minus sign here and no sign for ONe WDs?
        white_dwarfs['luminosity'] = -white_dwarfs.apply(
                estimate_by_magnitudes,
                axis=1,
                metallicity_grid=metallicities[spectral_type],
                cooling_sequences=cooling_sequences[spectral_type],
                pre_wd_lifetimes=pre_wd_lifetimes[spectral_type],
                interest_parameter='luminosity')
        white_dwarfs['effective_temperature'] = white_dwarfs.apply(
                estimate_by_magnitudes,
                axis=1,
                metallicity_grid=metallicities[spectral_type],
                cooling_sequences=cooling_sequences[spectral_type],
                pre_wd_lifetimes=pre_wd_lifetimes[spectral_type],
                interest_parameter='effective_temperature')

        for color in colors:
            white_dwarfs[color] = white_dwarfs.apply(
                    estimate_color,
                    axis=1,
                    color_table=color_tables[spectral_type])

    oxygen_neon_white_dwarfs['spectral_type'] = SpectralType.ONe
    parameters = ['luminosity',
                  'u_ubvri_absolute',
                  'b_ubvri_absolute',
                  'v_ubvri_absolute',
                  'r_ubvri_absolute',
                  'i_ubvri_absolute',
                  'j_ubvri_absolute',
                  'effective_temperature']

    for parameter in parameters:
        oxygen_neon_white_dwarfs.apply(estimate_oxygen_neon_parameters,
                                       axis=1,
                                       color_table=one_color_table,
                                       interest_parameter=parameter)

    return carbon_oxygen_white_dwarfs + oxygen_neon_white_dwarfs


def estimate_oxygen_neon_parameters(
        star: pd.Series,
        *,
        color_table: Dict[int, pd.DataFrame],
        interest_parameter: str) -> None:
    pre_wd_lifetime = 0.  # constant for all ONe white dwarfs

    mass = star['mass']
    cooling_time = star['cooling_time']
    mass_grid = np.array([key / 1e5
                          for key in color_table.keys()])
    int_mass_grid = list(color_table.keys())

    if mass < mass_grid[0] or mass >= mass_grid[-1]:
        estimate_interest_value = extrapolate_interest_value
    else:
        estimate_interest_value = interpolate_interest_value

    lesser_mass_index = calculate_index(mass,
                                        grid=mass_grid)
    lesser_int_mass = int_mass_grid[lesser_mass_index]
    lesser_mass_df = color_table[lesser_int_mass]
    greater_int_mass = int_mass_grid[lesser_mass_index + 1]
    greater_mass_df = color_table[greater_int_mass]

    estimate = partial(
            estimate_interest_value,
            mass=mass,
            cooling_time=cooling_time,
            greater_mass_cooling_time_grid=greater_mass_df['cooling_time'],
            greater_mass_pre_wd_lifetime=pre_wd_lifetime,
            lesser_mass_cooling_time_grid=greater_mass_df['cooling_time'],
            lesser_mass_pre_wd_lifetime=pre_wd_lifetime,
            min_mass=mass_grid[lesser_mass_index],
            max_mass=mass_grid[lesser_mass_index + 1])

    return estimate(greater_mass_interest_parameter_grid=greater_mass_df[
                        interest_parameter],
                    lesser_mass_interest_parameter_grid=lesser_mass_df[
                        interest_parameter])


def estimate_color(star: pd.Series,
                   *,
                   color_table: Dict[int, pd.DataFrame],
                   color: str) -> float:
    mass = star['mass']
    star_luminosity = star['luminosity']
    mass_grid = np.array([key / 1e5
                          for key in color_table.keys()])
    int_mass_grid = list(color_table.keys())

    lesser_mass_index = calculate_index(mass,
                                        grid=mass_grid)
    greater_mass_index = lesser_mass_index + 1
    lesser_int_mass = int_mass_grid[lesser_mass_index]
    lesser_mass_df = color_table[lesser_int_mass]
    greater_int_mass = int_mass_grid[greater_mass_index]
    greater_mass_df = color_table[greater_int_mass]

    min_luminosity_grid = lesser_mass_df['luminosity']
    max_luminosity_grid = greater_mass_df['luminosity']

    row_index = calculate_index(star_luminosity,
                                grid=min_luminosity_grid)
    next_row_index = calculate_index(
            star_luminosity,
            grid=max_luminosity_grid)

    if (star_luminosity > min_luminosity_grid[0] or
            star_luminosity > max_luminosity_grid[0]):
        min_mass = mass_grid[lesser_mass_index]
        max_mass = mass_grid[greater_mass_index]
    elif (star_luminosity < min_luminosity_grid[-1] or
            star_luminosity < max_luminosity_grid[-1]):
        min_mass = mass_grid[lesser_mass_index]
        max_mass = mass_grid[greater_mass_index]
    else:
        # TODO: check these indexes, they look suspicious
        min_mass = mass_grid[0]
        max_mass = mass_grid[1]

    min_magnitude_grid = lesser_mass_df[color]
    max_magnitude_grid = greater_mass_df[color]

    min_magnitude = estimate_at(
            star_luminosity,
            x=(min_luminosity_grid[row_index],
               min_luminosity_grid[row_index + 1]),
            y=(min_magnitude_grid[row_index],
               min_magnitude_grid[row_index + 1]))
    max_magnitude = estimate_at(
            star_luminosity,
            x=(max_luminosity_grid[next_row_index],
               max_luminosity_grid[next_row_index + 1]),
            y=(max_magnitude_grid[next_row_index],
               max_magnitude_grid[next_row_index + 1]))

    return estimate_at(mass,
                       x=(min_mass, max_mass),
                       y=(min_magnitude, max_magnitude))


def generate_spectral_types(*,
                            db_to_da_fraction: float,
                            size: int) -> np.ndarray:
    spectral_types = np.empty(size)

    randoms = np.random.rand(size)
    db_mask = randoms < db_to_da_fraction

    spectral_types[db_mask] = SpectralType.DB
    spectral_types[~db_mask] = SpectralType.DA

    return spectral_types


def estimate_by_magnitudes(
        star: pd.Series,
        *,
        metallicity_grid: List[float],
        cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        pre_wd_lifetimes: Dict[int, np.ndarray],
        interest_parameter: str) -> float:
    metallicity = star['metallicity']
    mass = star['mass']
    cooling_time = star['cooling_time']

    min_metallicity_index = get_min_metallicity_index(
            metallicity=metallicity,
            grid_metallicities=metallicity_grid)
    min_metallicity = metallicity_grid[min_metallicity_index]
    max_metallicity = metallicity_grid[min_metallicity_index + 1]

    int_min_metallicity = int(min_metallicity * 1e3)
    int_max_metallicity = int(max_metallicity * 1e3)

    min_metallicity_grids = cooling_sequences[int_min_metallicity]
    max_metallicity_grids = cooling_sequences[int_max_metallicity]

    min_metallicity_pre_wd_lifetimes = pre_wd_lifetimes[int_min_metallicity]
    max_metallicity_pre_wd_lifetimes = pre_wd_lifetimes[int_max_metallicity]

    estimate = partial(estimate_edge_case,
                       mass=mass,
                       cooling_time=cooling_time,
                       interest_sequence_str=interest_parameter)

    min_interest_parameter = estimate(
            cooling_sequences=min_metallicity_grids,
            pre_wd_lifetimes=min_metallicity_pre_wd_lifetimes)
    max_interest_parameter = estimate(
            cooling_sequences=max_metallicity_grids,
            pre_wd_lifetimes=max_metallicity_pre_wd_lifetimes)

    return estimate_at(metallicity,
                       x=(min_metallicity, max_metallicity),
                       y=(min_interest_parameter, max_interest_parameter))


# TODO: this looks like a case of ONe stars
def estimate_edge_case(*,
                       cooling_sequences: Dict[int, pd.DataFrame],
                       mass: float,
                       cooling_time: float,
                       interest_parameter: str,
                       pre_wd_lifetimes: np.ndarray) -> float:
    int_mass_grid = sorted(list(cooling_sequences.keys()))
    mass_grid = np.array([int_mass / 1e5
                          for int_mass in int_mass_grid])
    mass_index = calculate_index(mass,
                                 grid=mass_grid)
    int_lesser_mass_index = int_mass_grid[mass_index]
    int_greater_mass_index = int_mass_grid[mass_index + 1]
    lesser_mass_df = cooling_sequences[int_lesser_mass_index]
    greater_mass_df = cooling_sequences[int_greater_mass_index]
    lesser_mass_pre_wd_lifetime = pre_wd_lifetimes[int_lesser_mass_index]
    greater_mass_pre_wd_lifetime = pre_wd_lifetimes[int_greater_mass_index]

    if mass < mass_grid[0] or mass >= mass_grid[-1]:
        estimate_interest_value = extrapolate_interest_value
    else:
        estimate_interest_value = interpolate_interest_value

    return estimate_interest_value(
            mass=mass,
            cooling_time=cooling_time,
            greater_mass_cooling_time_grid=greater_mass_df['cooling_time'],
            greater_mass_interest_parameter_grid=greater_mass_df[
                interest_parameter],
            greater_mass_pre_wd_lifetime=greater_mass_pre_wd_lifetime,
            lesser_mass_cooling_time_grid=lesser_mass_df['cooling_time'],
            lesser_mass_interest_parameter_grid=lesser_mass_df[
                interest_parameter],
            lesser_mass_pre_wd_lifetime=lesser_mass_pre_wd_lifetime,
            min_mass=mass_grid[mass_index],
            max_mass=mass_grid[mass_index + 1])


def extrapolating_by_grid(cooling_time: float,
                          *,
                          cooling_time_grid: np.ndarray) -> bool:
    if (cooling_time < cooling_time_grid[0] or
            cooling_time >= cooling_time_grid[-1]):
        return True
    else:
        return False


def get_min_metallicity_index(*,
                              metallicity: float,
                              grid_metallicities: List[float]) -> int:
    if (metallicity < grid_metallicities[0] or
            metallicity > grid_metallicities[-1]):
        raise ValueError(f'There is no support for metallicities '
                         f'lying out of the range of {grid_metallicities}')
    metallicity = np.array([metallicity])
    left_index = np.searchsorted(grid_metallicities, metallicity) - 1.
    return np.asscalar(left_index)


def calculate_index(value: float,
                    *,
                    grid: np.ndarray) -> int:
    if value <= grid[0]:
        return 0
    elif value > grid[-1]:
        # Index of element before the last one
        return -2
    cooling_time = np.array([value])
    left_index = np.searchsorted(grid, cooling_time) - 1
    return np.asscalar(left_index)


def estimate_at(x_0: float,
                *,
                x: Tuple[float, float],
                y: Tuple[float, float]) -> float:
    spline = linear_estimation(x=x,
                               y=y)
    return spline(x_0)


def interpolate_interest_value(
        *,
        mass: float,
        cooling_time: float,
        greater_mass_cooling_time_grid: np.ndarray,
        greater_mass_interest_parameter_grid: np.ndarray,
        greater_mass_pre_wd_lifetime: float,
        lesser_mass_cooling_time_grid: np.ndarray,
        lesser_mass_interest_parameter_grid: np.ndarray,
        lesser_mass_pre_wd_lifetime: float,
        min_mass: float,
        max_mass: float) -> float:
    min_extrapolated_interest_value = partial(
            estimated_interest_value,
            cooling_time=cooling_time,
            cooling_time_grid=(lesser_mass_cooling_time_grid
                               + lesser_mass_pre_wd_lifetime),
            interest_sequence_grid=lesser_mass_interest_parameter_grid)
    max_extrapolated_interest_value = partial(
            estimated_interest_value,
            cooling_time=cooling_time,
            cooling_time_grid=(greater_mass_cooling_time_grid
                               + greater_mass_pre_wd_lifetime),
            interest_sequence_grid=greater_mass_interest_parameter_grid)

    extrapolating_by_min_cooling_time_grid = extrapolating_by_grid(
            cooling_time,
            cooling_time_grid=lesser_mass_cooling_time_grid)
    min_row_index = calculate_index(cooling_time,
                                    grid=lesser_mass_cooling_time_grid)

    if extrapolating_by_min_cooling_time_grid:
        x_1 = min_extrapolated_interest_value(row_index=min_row_index)
    else:
        y_1 = lesser_mass_cooling_time_grid[min_row_index]
        y_2 = lesser_mass_cooling_time_grid[min_row_index + 1]
        x_1 = lesser_mass_interest_parameter_grid[min_row_index]
        x_2 = lesser_mass_interest_parameter_grid[min_row_index + 1]

    extrapolating_by_max_cooling_time_grid = extrapolating_by_grid(
            cooling_time,
            cooling_time_grid=greater_mass_cooling_time_grid)
    max_row_index = calculate_index(cooling_time,
                                    grid=greater_mass_cooling_time_grid)
    if extrapolating_by_max_cooling_time_grid:
        x_3 = max_extrapolated_interest_value(row_index=max_row_index)
    else:
        y_3 = greater_mass_cooling_time_grid[max_row_index]
        y_4 = greater_mass_cooling_time_grid[max_row_index + 1]
        x_3 = greater_mass_interest_parameter_grid[max_row_index]
        x_4 = greater_mass_interest_parameter_grid[max_row_index + 1]

    if (not extrapolating_by_min_cooling_time_grid and
            not extrapolating_by_max_cooling_time_grid):
        ym_1 = estimate_at(mass,
                           x=(min_mass, max_mass),
                           y=(y_1, y_3))
        ym_2 = estimate_at(mass,
                           x=(min_mass, max_mass),
                           y=(y_2, y_4))
        xm_1 = estimate_at(mass,
                           x=(min_mass, max_mass),
                           y=(x_1, x_3))
        xm_2 = estimate_at(mass,
                           x=(min_mass, max_mass),
                           y=(x_2, x_4))

        return estimate_at(cooling_time,
                           x=(ym_1, ym_2),
                           y=(xm_1, xm_2))

    if (not extrapolating_by_min_cooling_time_grid
            and extrapolating_by_max_cooling_time_grid):
        xm_1 = estimate_at(cooling_time,
                           x=(y_1, y_2),
                           y=(x_1, x_2))
        return estimate_at(mass,
                           x=(min_mass, max_mass),
                           y=(xm_1, x_3))

    if (extrapolating_by_min_cooling_time_grid and
            not extrapolating_by_max_cooling_time_grid):
        xm_2 = estimate_at(cooling_time,
                           x=(y_3, y_4),
                           y=(x_3, x_4))
        return estimate_at(mass,
                           x=(min_mass, max_mass),
                           y=(x_1, xm_2))

    return estimate_at(mass,
                       x=(min_mass, max_mass),
                       y=(x_1, x_3))


def extrapolate_interest_value(
        *,
        mass: float,
        cooling_time: float,
        greater_mass_cooling_time_grid: np.ndarray,
        greater_mass_interest_parameter_grid: np.ndarray,
        greater_mass_pre_wd_lifetime: float,
        lesser_mass_cooling_time_grid: np.ndarray,
        lesser_mass_interest_parameter_grid: np.ndarray,
        lesser_mass_pre_wd_lifetime: float,
        min_mass: float,
        max_mass: float) -> float:
    interest_value = partial(get_interest_value,
                             cooling_time=cooling_time)
    min_interest_value = interest_value(
            cooling_time_grid=lesser_mass_cooling_time_grid,
            interest_sequence_grid=lesser_mass_interest_parameter_grid,
            pre_wd_lifetime=lesser_mass_pre_wd_lifetime)
    max_interest_value = interest_value(
            cooling_time_grid=greater_mass_cooling_time_grid,
            interest_sequence_grid=greater_mass_interest_parameter_grid,
            pre_wd_lifetime=greater_mass_pre_wd_lifetime)
    return estimate_at(mass,
                       x=(min_mass, max_mass),
                       y=(min_interest_value, max_interest_value))


def get_interest_value(*,
                       cooling_time: float,
                       cooling_time_grid: np.ndarray,
                       pre_wd_lifetime: float,
                       interest_parameter_grid: np.ndarray) -> float:
    row_index = calculate_index(cooling_time,
                                grid=cooling_time_grid)

    if (cooling_time < cooling_time_grid[0] or
            cooling_time > cooling_time_grid[-1]):
        cooling_time_grid = cooling_time_grid + pre_wd_lifetime

    return estimated_interest_value(
            row_index=row_index,
            cooling_time=cooling_time,
            cooling_time_grid=cooling_time_grid,
            interest_parameter_grid=interest_parameter_grid)


def estimated_interest_value(*,
                             cooling_time: float,
                             cooling_time_grid: np.ndarray,
                             interest_parameter_grid: np.ndarray,
                             row_index: int) -> float:
    return estimate_at(cooling_time,
                       x=(cooling_time_grid[row_index],
                          cooling_time_grid[row_index + 1]),
                       y=(interest_parameter_grid[row_index],
                          interest_parameter_grid[row_index + 1]))

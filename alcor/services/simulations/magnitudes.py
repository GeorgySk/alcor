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
        one_color_table: Dict[int, pd.DataFrame]) -> pd.DataFrame:
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
                            == SpectralType.DA.value)
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
        if white_dwarfs.shape[0] == 0:
            continue
        # TODO: why is there minus sign here and no sign for ONe WDs?
        white_dwarfs.loc[:, 'luminosity'] = -white_dwarfs.apply(
                estimate_by_metallicities,
                axis=1,
                metallicity_grid=metallicities[spectral_type],
                cooling_sequences=cooling_sequences[spectral_type],
                interest_parameter='luminosity')
        white_dwarfs.loc[:, 'effective_temperature'] = white_dwarfs.apply(
                estimate_by_metallicities,
                axis=1,
                metallicity_grid=metallicities[spectral_type],
                cooling_sequences=cooling_sequences[spectral_type],
                interest_parameter='effective_temperature')

        for color in colors:
            if white_dwarfs.shape[0] == 0:
                continue
            white_dwarfs.loc[:, color] = white_dwarfs.apply(
                    estimate_color,
                    axis=1,
                    color_table=color_tables[spectral_type],
                    color=color)
    oxygen_neon_white_dwarfs['spectral_type'] = SpectralType.ONe.value
    parameters = ['luminosity',
                  'u_ubvri_absolute',
                  'b_ubvri_absolute',
                  'v_ubvri_absolute',
                  'r_ubvri_absolute',
                  'i_ubvri_absolute',
                  'j_ubvri_absolute',
                  'effective_temperature']

    for parameter in parameters:
        oxygen_neon_white_dwarfs.apply(estimate_by_mass,
                                       axis=1,
                                       color_table=one_color_table,
                                       interest_parameter=parameter)

    return pd.concat([da_white_dwarfs,
                      db_white_dwarfs,
                      oxygen_neon_white_dwarfs])


def estimate_by_metallicities(
        star: pd.Series,
        *,
        metallicity_grid: List[float],
        cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        interest_parameter: str) -> float:
    metallicity = star['metallicity']

    min_metallicity_index = get_min_metallicity_index(
            metallicity=metallicity,
            grid_metallicities=metallicity_grid)
    min_metallicity = metallicity_grid[min_metallicity_index]
    max_metallicity = metallicity_grid[min_metallicity_index + 1]

    int_min_metallicity = int(min_metallicity * 1e3)
    int_max_metallicity = int(max_metallicity * 1e3)

    min_metallicity_grids = cooling_sequences[int_min_metallicity]
    max_metallicity_grids = cooling_sequences[int_max_metallicity]

    estimate = partial(estimate_by_mass,
                       star=star,
                       interest_parameter=interest_parameter)

    min_interest_parameter = estimate(tracks=min_metallicity_grids)
    max_interest_parameter = estimate(tracks=max_metallicity_grids)

    return estimate_at(metallicity,
                       x=(min_metallicity, max_metallicity),
                       y=(min_interest_parameter, max_interest_parameter))


def estimate_by_mass(
        star: pd.Series,
        *,
        tracks: Dict[int, pd.DataFrame],
        interest_parameter: str) -> float:
    mass = star['mass']
    cooling_time = star['cooling_time']

    int_mass_grid = sorted(list(tracks.keys()))
    mass_grid = np.array([key / 1e5
                          for key in int_mass_grid])

    lesser_mass_index = calculate_index(mass,
                                        grid=mass_grid)
    lesser_int_mass = int_mass_grid[lesser_mass_index]
    lesser_mass_df = tracks[lesser_int_mass]
    greater_mass_index = lesser_mass_index + 1
    greater_int_mass = int_mass_grid[greater_mass_index]
    greater_mass_df = tracks[greater_int_mass]

    if mass < mass_grid[0] or mass >= mass_grid[-1]:
        estimate_interest_value = extrapolate_interest_value
    else:
        estimate_interest_value = interpolate_interest_value

    min_row_index = calculate_index(cooling_time,
                                    grid=lesser_mass_df['cooling_time'].values)
    max_row_index = calculate_index(
            cooling_time,
            grid=greater_mass_df['cooling_time'].values)

    return estimate_interest_value(
            mass=mass,
            cooling_time=cooling_time,
            greater_mass_cooling_time_grid=(
                greater_mass_df['cooling_time'].values),
            greater_mass_interest_parameter_grid=greater_mass_df[
                interest_parameter].values,
            lesser_mass_cooling_time_grid=(
                lesser_mass_df['cooling_time'].values),
            lesser_mass_interest_parameter_grid=lesser_mass_df[
                interest_parameter].values,
            min_mass=mass_grid[lesser_mass_index],
            max_mass=mass_grid[greater_mass_index],
            min_row_index=min_row_index,
            max_row_index=max_row_index)


def estimate_color(star: pd.Series,
                   *,
                   color_table: Dict[int, pd.DataFrame],
                   color: str) -> float:
    mass = star['mass']
    star_luminosity = star['luminosity']

    int_mass_grid = sorted(list(color_table.keys()))
    mass_grid = np.array([key / 1e5
                          for key in int_mass_grid])

    lesser_mass_index = calculate_index(mass,
                                        grid=mass_grid)
    lesser_int_mass = int_mass_grid[lesser_mass_index]
    lesser_mass_df = color_table[lesser_int_mass]
    greater_mass_index = lesser_mass_index + 1
    greater_int_mass = int_mass_grid[greater_mass_index]
    greater_mass_df = color_table[greater_int_mass]

    min_luminosity_grid = lesser_mass_df['luminosity'].values
    max_luminosity_grid = greater_mass_df['luminosity'].values

    row_index = calculate_index(star_luminosity,
                                grid=min_luminosity_grid)
    next_row_index = calculate_index(
            star_luminosity,
            grid=max_luminosity_grid)

    min_mass = mass_grid[lesser_mass_index]
    max_mass = mass_grid[greater_mass_index]

    min_magnitude_grid = lesser_mass_df[color].values
    max_magnitude_grid = greater_mass_df[color].values

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


def interpolate_interest_value(
        *,
        mass: float,
        cooling_time: float,
        greater_mass_cooling_time_grid: np.ndarray,
        greater_mass_interest_parameter_grid: np.ndarray,
        lesser_mass_cooling_time_grid: np.ndarray,
        lesser_mass_interest_parameter_grid: np.ndarray,
        min_mass: float,
        max_mass: float,
        min_row_index: int,
        max_row_index: int) -> float:
    extrapolating_by_min_cooling_time_grid = extrapolating_by_grid(
            cooling_time,
            cooling_time_grid=lesser_mass_cooling_time_grid)
    if extrapolating_by_min_cooling_time_grid:
        x_1 = estimated_interest_value(
                row_index=min_row_index,
                cooling_time=cooling_time,
                cooling_time_grid=lesser_mass_cooling_time_grid,
                interest_parameter_grid=lesser_mass_interest_parameter_grid)
    else:
        y_1 = lesser_mass_cooling_time_grid[min_row_index]
        y_2 = lesser_mass_cooling_time_grid[min_row_index + 1]
        x_1 = lesser_mass_interest_parameter_grid[min_row_index]
        x_2 = lesser_mass_interest_parameter_grid[min_row_index + 1]

    extrapolating_by_max_cooling_time_grid = extrapolating_by_grid(
            cooling_time,
            cooling_time_grid=greater_mass_cooling_time_grid)
    if extrapolating_by_max_cooling_time_grid:
        x_3 = estimated_interest_value(
                row_index=max_row_index,
                cooling_time=cooling_time,
                cooling_time_grid=greater_mass_cooling_time_grid,
                interest_parameter_grid=greater_mass_interest_parameter_grid)
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
        lesser_mass_cooling_time_grid: np.ndarray,
        lesser_mass_interest_parameter_grid: np.ndarray,
        min_mass: float,
        max_mass: float,
        min_row_index: int,
        max_row_index: int) -> float:
    min_interest_value = estimated_interest_value(
            row_index=min_row_index,
            cooling_time=cooling_time,
            cooling_time_grid=lesser_mass_cooling_time_grid,
            interest_parameter_grid=lesser_mass_interest_parameter_grid)
    max_interest_value = estimated_interest_value(
            row_index=max_row_index,
            cooling_time=cooling_time,
            cooling_time_grid=greater_mass_cooling_time_grid,
            interest_parameter_grid=greater_mass_interest_parameter_grid)

    return estimate_at(mass,
                       x=(min_mass, max_mass),
                       y=(min_interest_value, max_interest_value))


def generate_spectral_types(*,
                            db_to_da_fraction: float,
                            size: int) -> np.ndarray:
    spectral_types = np.empty(size,
                              dtype='<U3')

    randoms = np.random.rand(size)
    db_mask = randoms < db_to_da_fraction

    # TODO: use SpectralType.DB.value
    spectral_types[db_mask] = SpectralType.DB.value
    spectral_types[~db_mask] = SpectralType.DA.value

    return spectral_types


def get_min_metallicity_index(metallicity: float,
                              *,
                              grid_metallicities: List[float]) -> int:
    if (metallicity < grid_metallicities[0] or
            metallicity > grid_metallicities[-1]):
        raise ValueError('There is no support for metallicities '
                         'lying out of the range of {grid_metallicities}'
                         .format(grid_metallicities=grid_metallicities))
    metallicity = np.array([metallicity])
    left_index = np.searchsorted(grid_metallicities, metallicity) - 1.
    left_index = left_index.astype(int)
    left_index = np.asscalar(left_index)

    if metallicity == grid_metallicities[left_index + 1]:
        return left_index + 1

    return left_index


def extrapolating_by_grid(cooling_time: float,
                          *,
                          cooling_time_grid: np.ndarray) -> bool:
    if (cooling_time < cooling_time_grid[0] or
            cooling_time >= cooling_time_grid[-1]):
        return True
    else:
        return False


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
                x: Tuple[float, ...],
                y: Tuple[float, ...]) -> float:
    spline = linear_estimation(x=x,
                               y=y)
    return np.asscalar(spline(x_0))

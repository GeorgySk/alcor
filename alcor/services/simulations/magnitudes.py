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


def assign_magnitudes(stars: pd.DataFrame,
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

        da_db_interpolation(
                star,
                cooling_sequences=cooling_sequences[spectral_type],
                color_table=colo_tables[spectral_type],
                metallicities=metallicities[spectral_type])

    oxygen_neon_white_dwarfs['cooling_time'] = 9. + log10(
            oxygen_neon_white_dwarfs['cooling_time'])
    oxygen_neon_white_dwarfs['spectral_type'] = SpectralType.ONe

    for _, star in oxygen_neon_white_dwarfs.iterrows():
        one_interpolation(star,
                          color_table=one_color_table)

    return carbon_oxygen_white_dwarfs + oxygen_neon_white_dwarfs


def generate_spectral_type(db_to_da_fraction: float) -> SpectralType:
    if np.random.rand() < db_to_da_fraction:
        return SpectralType.DB
    return SpectralType.DA


def one_interpolation(star: pd.Series,
                      *,
                      color_table: Dict[str, np.ndarray],
                      one_model: bool = True,
                      by_logarithm: bool = False) -> None:
    star_mass = star['mass']
    star_cooling_time = star['cooling_time']
    mass_grid = color_table['mass']
    cooling_time_grid = color_table['cooling_time']
    pre_wd_lifetime_grid = color_table['pre_wd_lifetime_grid']
    rows_counts = color_table['rows_counts']

    extrapolate = partial(extrapolate_interest_value,
                          star_mass=star_mass,
                          star_cooling_time=star_cooling_time,
                          cooling_time_grid=cooling_time_grid,
                          pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                          rows_counts=rows_counts,
                          by_logarithm=by_logarithm)

    if star_mass < mass_grid[0]:
        do_estimation = partial(extrapolate,
                                min_mass_index=0,
                                min_mass=mass_grid[0],
                                max_mass=mass_grid[1])
    elif star_mass >= mass_grid[-1]:
        min_mass_index = mass_grid.size() - 1
        do_estimation = partial(extrapolate,
                                min_mass_index=min_mass_index,
                                min_mass=mass_grid[min_mass_index],
                                max_mass=mass_grid[min_mass_index + 1])
    else:
        min_mass_index = get_mass_index(star_mass=star_mass,
                                        mass_grid=mass_grid)
        min_mass = mass_grid[min_mass_index]
        max_mass = mass_grid[min_mass_index + 1]
        do_estimation = partial(interpolate_by_mass,
                                star=star,
                                cooling_or_color_sequence=color_table,
                                min_mass=min_mass,
                                max_mass=max_mass,
                                min_mass_index=min_mass_index,
                                by_logarithm=by_logarithm,
                                one_model=one_model)

    star['luminosity'] = do_estimation(interest_sequence_grid='luminosity')
    v_ubvri_absolute = do_estimation(
            interest_sequence_grid='v_ubvri_absolute')
    bv_ubvri = do_estimation(interest_sequence_grid='bv_ubvri')
    vi_ubvri = do_estimation(interest_sequence_grid='vi_ubvri')
    vr_ubvri = do_estimation(interest_sequence_grid='vr_ubvri')
    uv_ubvri = do_estimation(interest_sequence_grid='uv_ubvri')
    log_effective_temperature = do_estimation(
            interest_sequence_grid='log_effective_temperature')

    star['effective_temperature'] = 10. ** log_effective_temperature

    star['u_ubvri_absolute'] = uv_ubvri + v_ubvri_absolute
    star['b_ubvri_absolute'] = bv_ubvri + v_ubvri_absolute
    star['r_ubvri_absolute'] = v_ubvri_absolute - vr_ubvri
    star['i_ubvri_absolute'] = v_ubvri_absolute - vi_ubvri
    star['v_ubvri_absolute'] = v_ubvri_absolute


def get_min_metallicity_index(*,
                              star_metallicity: float,
                              grid_metallicities: List[float]) -> int:
    if (star_metallicity < grid_metallicities[0]
            or star_metallicity > grid_metallicities[-1]):
        raise ValueError(f'There is no support for metallicities '
                         f'lying out of the range of {grid_metallicities}')
    star_metallicity = np.array([star_metallicity])
    left_index = np.searchsorted(grid_metallicities, star_metallicity) - 1.
    return np.asscalar(left_index)


def estimate_edge_case(*,
                       star_mass: float,
                       star_cooling_time: float,
                       cooling_time_grid: np.ndarray,
                       mass_grid: np.ndarray,
                       pre_wd_lifetime_grid: np.ndarray,
                       rows_counts: np.ndarray,
                       interest_sequence_grid: np.ndarray,
                       by_logarithm: bool) -> float:
    if star_mass < mass_grid[0]:
        min_mass_index = 0
        do_estimation = extrapolate_interest_value
    elif star_mass >= mass_grid[-1]:
        min_mass_index = mass_grid.size() - 1
        do_estimation = extrapolate_interest_value
    else:
        min_mass_index = get_mass_index(star_mass=star_mass,
                                        mass_grid=mass_grid)
        do_estimation = interpolate_by_mass
    return do_estimation(
            star_mass=star_mass,
            star_cooling_time=star_cooling_time,
            cooling_time_grid=cooling_time_grid,
            pre_wd_lifetime_grid=pre_wd_lifetime_grid,
            rows_counts=rows_counts,
            min_mass_index=min_mass_index,
            interest_sequence_grid=interest_sequence_grid,
            min_mass=mass_grid[min_mass_index],
            max_mass=mass_grid[min_mass_index + 1],
            by_logarithm=by_logarithm)


def da_db_interpolation(star: pd.Series,
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

    estimate_min_case = partial(
            estimate_edge_case,
            star_mass=star_mass,
            star_cooling_time=star_cooling_time,
            cooling_time_grid=min_metallicity_grids['cooling_time'],
            mass_grid=min_metallicity_grids['mass'],
            pre_wd_lifetime_grid=min_metallicity_grids['pre_wd_lifetime'],
            rows_counts=min_metallicity_grids['rows_counts'])
    estimate_max_case = partial(
            estimate_edge_case,
            star_mass=star_mass,
            star_cooling_time=star_cooling_time,
            cooling_time_grid=max_metallicity_grids['cooling_time'],
            mass_grid=max_metallicity_grids['mass'],
            pre_wd_lifetime_grid=max_metallicity_grids['pre_wd_lifetime'],
            rows_counts=max_metallicity_grids['rows_counts'])

    min_luminosity = estimate_min_case(
            interest_sequence_grid=min_metallicity_grids['luminosity'],
            by_logarithm=False)
    max_luminosity = estimate_max_case(
            interest_sequence_grid=max_metallicity_grids['luminosity'],
            by_logarithm=False)
    min_effective_temperature = estimate_min_case(
            interest_sequence_grid=min_metallicity_grids[
                'effective_temperature'],
            by_logarithm=True)
    max_effective_temperature = estimate_max_case(
            interest_sequence_grid=max_metallicity_grids[
                'effective_temperature'],
            by_logarithm=True)

    star['luminosity'] = -estimate_at(star_metallicity,
                                      x=(min_metallicity, max_metallicity),
                                      y=(min_luminosity, max_luminosity))

    star['effective_temperature'] = estimate_at(
            star_metallicity,
            x=(min_metallicity, max_metallicity),
            y=(min_effective_temperature, max_effective_temperature))

    star = star_with_colors(star,
                            color_table=color_table)

    return star


def star_with_colors(star: pd.Series,
                     *,
                     color_table: Dict[str, np.ndarray]) -> pd.Series:
    star_mass = star['mass']
    star_luminosity = star['luminosity']
    luminosity_grid = color_table['luminosity']
    mass_grid = color_table['mass_grid']
    rows_counts = color_table['rows_counts']

    min_mass_index = find_mass_index(star_mass=star_mass,
                                     mass_grid=mass_grid)
    rows_count = rows_counts[min_mass_index]
    next_rows_count = rows_counts[min_mass_index + 1]

    colors = ['u_ubvri_absolute',
              'b_ubvri_absolute',
              'v_ubvri_absolute',
              'r_ubvri_absolute',
              'i_ubvri_absolute']

    if (star_luminosity > luminosity_grid[min_mass_index, 0]
            or star_luminosity > luminosity_grid[min_mass_index + 1, 0]):
        min_mass = mass_grid[min_mass_index]
        max_mass = mass_grid[min_mass_index + 1]
        row_index = 0
        next_row_index = 0
    elif (star_luminosity < luminosity_grid[min_mass_index, rows_count]
          or star_luminosity < luminosity_grid[min_mass_index + 1,
                                               next_rows_count]):
        min_mass = mass_grid[min_mass_index]
        max_mass = mass_grid[min_mass_index + 1]
        row_index = rows_count
        next_row_index = next_rows_count
    else:
        min_mass = mass_grid[0]
        max_mass = mass_grid[1]
        find_row_index = partial(find_index,
                                 luminosity=star_luminosity)
        row_index = find_row_index(
                rows_count=rows_count,
                luminosity_grid_for_specific_mass=luminosity_grid[
                                                  min_mass_index, :])
        next_row_index = find_row_index(
                rows_count=next_rows_count,
                luminosity_grid_for_specific_mass=luminosity_grid[
                                                  min_mass_index + 1, :])

    for color in colors:
        magnitude_grid = color_table[color]
        min_magnitude = estimate_at(
                star_luminosity,
                x=(luminosity_grid[min_mass_index, row_index],
                   luminosity_grid[min_mass_index, row_index + 1]),
                y=(magnitude_grid[min_mass_index, row_index],
                   magnitude_grid[min_mass_index, row_index + 1]))
        max_magnitude = estimate_at(
                star_luminosity,
                x=(luminosity_grid[min_mass_index + 1, next_row_index],
                   luminosity_grid[min_mass_index + 1, next_row_index + 1]),
                y=(magnitude_grid[min_mass_index + 1, next_row_index],
                   magnitude_grid[min_mass_index + 1, next_row_index + 1]))
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
                               rows_counts: np.ndarray,
                               min_mass: float,
                               max_mass: float,
                               by_logarithm: bool,
                               one_model: bool = False) -> float:
    interest_value = partial(get_interest_value,
                             star_cooling_time=star_cooling_time,
                             cooling_time_grid=cooling_time_grid,
                             pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                             interest_sequence_grid=interest_sequence_grid,
                             by_logarithm=by_logarithm,
                             one_model=one_model)

    min_interest_value = interest_value(
            mass_index=min_mass_index,
            rows_count=rows_counts[min_mass_index])
    max_interest_value = interest_value(
            mass_index=min_mass_index + 1,
            rows_count=rows_counts[min_mass_index + 1])
    return estimate_at(star_mass,
                       x=(min_mass, max_mass),
                       y=(min_interest_value, max_interest_value))


def get_mass_index(*,
                   star_mass: float,
                   mass_grid: np.ndarray) -> int:
    for row_index in range(mass_grid.size() - 1):
        if mass_grid[row_index] <= star_mass < mass_grid[row_index + 1]:
            return row_index


def interpolate_by_mass(*,
                        star_mass: float,
                        star_cooling_time: float,
                        min_mass: float,
                        max_mass: float,
                        min_mass_index: int,
                        cooling_time_grid: np.ndarray,
                        pre_wd_lifetime_grid: np.ndarray,
                        interest_sequence_grid: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    max_mass_index = min_mass_index + 1

    if one_model:
        extrapolated_interest_value = partial(
                one_white_dwarfs_estimated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                interest_sequence_grid=interest_sequence_grid)
    elif by_logarithm:
        extrapolated_interest_value = partial(
                extrapolated_interest_value_by_log,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                interest_sequence_grid=interest_sequence_grid)
    else:
        extrapolated_interest_value = partial(
                get_extrapolated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                interest_sequence_grid=interest_sequence_grid,
                by_logarithm=by_logarithm,
                one_model=one_model)

    if star_cooling_time < cooling_time_grid[min_mass_index, 0]:
        x1 = extrapolated_interest_value(min_row_index=1,
                                         mass_index=min_mass_index)
        case_1 = 1
    elif star_cooling_time >= cooling_time_grid[min_mass_index,
                                                rows_counts[min_mass_index]]:
        rows_count = rows_counts[min_mass_index]
        x1 = extrapolated_interest_value(min_row_index=rows_count,
                                         mass_index=min_mass_index)
        case_1 = 1
    else:
        for row_index in range(rows_counts[min_mass_index] - 1):
            if (cooling_time_grid[min_mass_index, row_index]
                    <= star_cooling_time
                    <= cooling_time_grid[min_mass_index, row_index + 1]):
                y1 = cooling_time_grid[min_mass_index, row_index]
                y2 = cooling_time_grid[min_mass_index, row_index + 1]
                x1 = interest_sequence_grid[min_mass_index, row_index]
                x2 = interest_sequence_grid[min_mass_index, row_index + 1]
                case_1 = 0

    if star_cooling_time < cooling_time_grid[max_mass_index, 0]:
        x3 = extrapolated_interest_value(min_row_index=1,
                                         mass_index=max_mass_index)
        case_2 = 1
    elif star_cooling_time >= cooling_time_grid[max_mass_index,
                                                rows_counts[max_mass_index]]:
        rows_count = rows_counts[max_mass_index]
        x3 = extrapolated_interest_value(min_row_index=rows_count,
                                         mass_index=max_mass_index)
        case_2 = 1
    else:
        for row_index in range(rows_counts[max_mass_index] - 1):
            if (cooling_time_grid[max_mass_index, row_index]
                    <= star_cooling_time
                    <= cooling_time_grid[max_mass_index, row_index + 1]):
                y3 = cooling_time_grid[max_mass_index, row_index]
                y4 = cooling_time_grid[max_mass_index, row_index + 1]
                x3 = interest_sequence_grid[max_mass_index, row_index]
                x4 = interest_sequence_grid[max_mass_index, row_index + 1]
                case_2 = 0

    if case_1 == 0 and case_2 == 0:
        ym1 = estimate_at(star_mass,
                          x=(min_mass, max_mass),
                          y=(y1, y3))
        ym2 = estimate_at(star_mass,
                          x=(min_mass, max_mass),
                          y=(y2, y4))
        xm1 = estimate_at(star_mass,
                          x=(min_mass, max_mass),
                          y=(x1, x3))
        xm2 = estimate_at(star_mass,
                          x=(min_mass, max_mass),
                          y=(x2, x4))

        return estimate_at(star_cooling_time,
                           x=(ym1, ym2),
                           y=(xm1, xm2))

    if case_1 == 0 and case_2 == 1:
        xm1 = estimate_at(star_cooling_time,
                          x=(y1, y2),
                          y=(x1, x2))
        return estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(xm1, x3))

    if case_1 == 1 and case_2 == 0:
        xm2 = estimate_at(star_cooling_time,
                          x=(y3, y4),
                          y=(x3, x4))
        return estimate_at(star_mass,
                           x=(min_mass, max_mass),
                           y=(x1, xm2))

    return estimate_at(star_mass,
                       x=(min_mass, max_mass),
                       y=(x1, x3))


def get_interest_value(*,
                       star_cooling_time: float,
                       cooling_time_grid: np.ndarray,
                       pre_wd_lifetime_grid: np.ndarray,
                       interest_sequence_grid: np.ndarray,
                       rows_count: int,
                       mass_index: int,
                       by_logarithm: bool,
                       one_model: bool = False) -> float:
    if one_model:
        extrapolated_interest_value = partial(
                one_white_dwarfs_estimated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                interest_sequence_grid=interest_sequence_grid,
                mass_index=mass_index)
    elif by_logarithm:
        extrapolated_interest_value = partial(
                extrapolated_interest_value_by_log,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                interest_sequence_grid=interest_sequence_grid,
                mass_index=mass_index)
    else:
        extrapolated_interest_value = partial(
                get_extrapolated_interest_value,
                star_cooling_time=star_cooling_time,
                cooling_time_grid=cooling_time_grid,
                pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                interest_sequence_grid=interest_sequence_grid,
                mass_index=mass_index,
                by_logarithm=by_logarithm,
                one_model=one_model)

    if star_cooling_time < cooling_time_grid[mass_index, 0]:
        return extrapolated_interest_value(min_row_index=0)

    if star_cooling_time > cooling_time_grid[mass_index, rows_count]:
        return extrapolated_interest_value(min_row_index=rows_count - 1)

    for row_index in range(rows_count - 1):
        if (cooling_time_grid[mass_index, row_index] <= star_cooling_time
                <= cooling_time_grid[mass_index, row_index + 1]):
            return estimate_at(
                    star_cooling_time,
                    x=(cooling_time_grid[mass_index, row_index],
                       cooling_time_grid[mass_index, row_index + 1]),
                    y=(interest_sequence_grid[mass_index, row_index],
                       interest_sequence_grid[mass_index, row_index + 1]))


def one_white_dwarfs_estimated_interest_value(
        *,
        star_cooling_time: float,
        cooling_time_grid: np.ndarray,
        interest_sequence_grid: np.ndarray,
        mass_index: int,
        min_row_index: int) -> float:
    return estimate_at(
            star_cooling_time,
            x=(cooling_time_grid[mass_index, min_row_index],
               cooling_time_grid[mass_index, min_row_index + 1]),
            y=(interest_sequence_grid[mass_index, min_row_index],
               interest_sequence_grid[mass_index, min_row_index + 1]))


def extrapolated_interest_value_by_log(
        *,
        star_cooling_time: float,
        cooling_time_grid: np.ndarray,
        interest_sequence_grid: np.ndarray,
        pre_wd_lifetime_grid: np.ndarray,
        mass_index: int,
        min_row_index: int) -> float:
    return 10. ** estimate_at(
            star_cooling_time,
            x=(log10(cooling_time_grid[mass_index, min_row_index]
                     + pre_wd_lifetime_grid[mass_index]),
               log10(cooling_time_grid[mass_index, min_row_index + 1]
                     + pre_wd_lifetime_grid[mass_index])),
            y=(log10(interest_sequence_grid[mass_index, min_row_index]),
               log10(interest_sequence_grid[mass_index, min_row_index + 1])))


def get_extrapolated_interest_value(*,
                                    star_cooling_time: float,
                                    cooling_time_grid: np.ndarray,
                                    pre_wd_lifetime_grid: np.ndarray,
                                    interest_sequence_grid: np.ndarray,
                                    mass_index: int,
                                    min_row_index: int) -> float:
    return estimate_at(
            star_cooling_time,
            x=(log10(cooling_time_grid[mass_index, min_row_index]
                     + pre_wd_lifetime_grid[mass_index]),
               log10(cooling_time_grid[mass_index, min_row_index + 1]
                     + pre_wd_lifetime_grid[mass_index])),
            y=(interest_sequence_grid[mass_index, min_row_index],
               interest_sequence_grid[mass_index, min_row_index + 1]))


def find_mass_index(*,
                    star_mass: float,
                    mass_grid: np.ndarray) -> int:
    if star_mass <= mass_grid[0]:
        return 0
    elif star_mass > mass_grid[-1]:
        # Index of element before the last one
        return -2
    star_mass = np.array([star_mass])
    left_index = np.searchsorted(mass_grid, star_mass) - 1
    return np.asscalar(left_index)


def find_index(*,
               rows_count: int,
               star_luminosity: float,
               luminosity_grid_for_specific_mass: np.ndarray) -> int:
    for row_index in range(rows_count - 1):
        if (luminosity_grid_for_specific_mass[row_index + 1]
                <= star_luminosity
                <= luminosity_grid_for_specific_mass[row_index]):
            return row_index


def estimate_at(x0: float,
                *,
                x: Tuple[float, float],
                y: Tuple[float, float]) -> float:
    spline = linear_estimation(x=x,
                               y=y)
    return spline(x0)

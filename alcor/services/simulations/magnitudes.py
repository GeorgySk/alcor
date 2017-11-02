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
                # TODO. can they be taken from cool.seq. keys?
                metallicities=metallicities[spectral_type])

    oxygen_neon_white_dwarfs['cooling_time'] = 9. + log10(
            oxygen_neon_white_dwarfs['cooling_time'])

    for _, star in oxygen_neon_white_dwarfs.iterrows():
        star['spectral_type'] = SpectralType.ONe
        (luminosity,
         effective_temperature,
         u_ubvri_absolute,
         b_ubvri_absolute,
         v_ubvri_absolute,
         r_ubvri_absolute,
         i_ubvri_absolute) = one_interpolation(star=star,
                                               color_table=one_color_table)

        star['luminosity'] = -luminosity
        star['effective_temperature'] = effective_temperature
        star['u_ubvri_absolute'] = u_ubvri_absolute
        star['b_ubvri_absolute'] = b_ubvri_absolute
        star['v_ubvri_absolute'] = v_ubvri_absolute
        star['r_ubvri_absolute'] = r_ubvri_absolute
        star['i_ubvri_absolute'] = i_ubvri_absolute

    return carbon_oxygen_white_dwarfs + oxygen_neon_white_dwarfs


def generate_spectral_type(db_to_da_fraction: float) -> SpectralType:
    if np.random.rand() < db_to_da_fraction:
        return SpectralType.DB
    return SpectralType.DA


def one_interpolation(*,
                      star: pd.Series,
                      color_table: Dict[str, np.ndarray],
                      one_model: bool = True,
                      by_logarithm: bool = False) -> Tuple[float, ...]:
    do_interpolation = partial(interpolate,
                               star=star,
                               cooling_or_color_sequence=color_table,
                               by_logarithm=by_logarithm,
                               one_model=one_model)
    luminosity = do_interpolation(interest_sequence='luminosity')
    v_ubvri_absolute = do_interpolation(interest_sequence='v_ubvri_absolute')
    bv_ubvri = do_interpolation(interest_sequence='bv_ubvri')
    vi_ubvri = do_interpolation(interest_sequence='vi_ubvri')
    vr_ubvri = do_interpolation(interest_sequence='vr_ubvri')
    uv_ubvri = do_interpolation(interest_sequence='uv_ubvri')
    log_effective_temperature = do_interpolation(
            interest_sequence='log_effective_temperature')

    effective_temperature = 10. ** log_effective_temperature

    u_ubvri_absolute = uv_ubvri + v_ubvri_absolute
    b_ubvri_absolute = bv_ubvri + v_ubvri_absolute
    r_ubvri_absolute = v_ubvri_absolute - vr_ubvri
    i_ubvri_absolute = v_ubvri_absolute - vi_ubvri

    return (luminosity,
            effective_temperature,
            u_ubvri_absolute,
            b_ubvri_absolute,
            v_ubvri_absolute,
            r_ubvri_absolute,
            i_ubvri_absolute)


def get_min_metallicity_index(*,
                              star_metallicity: float,
                              grid_metallicities: List[float]) -> int:
    for metallicity_index in range(len(grid_metallicities) - 1):
        if (grid_metallicities[metallicity_index]
                <= star_metallicity
                < grid_metallicities[metallicity_index + 1]):
            return metallicity_index


def da_db_interpolation(star: pd.Series,
                        *,
                        cooling_sequences: Dict[int, Dict[str, np.ndarray]],
                        color_table: Dict[str, np.ndarray],
                        metallicities: List[float]) -> None:
    luminosity_grid = color_table['luminosity']
    mass_grid = color_table['mass_grid']
    rows_counts = color_table['rows_counts']
    star_metallicity = star['metallicity']
    star_luminosity = star['luminosity']
    star_mass = star['mass']
    star_cooling_time = star['cooling_time']

    min_metallicity_index = get_min_metallicity_index(
            star_metallicity=star_metallicity,
            grid_metallicities=metallicities)
    min_metallicity = metallicities[min_metallicity_index]
    max_metallicity = metallicities[min_metallicity_index + 1]

    (min_luminosity,
     max_luminosity,
     min_effective_temperature,
     max_effective_temperature) = (
        get_luminosity_effective_temperature_limits(
            star_mass=star_mass,
            star_cooling_time=star_cooling_time,
            cooling_sequences=cooling_sequences,
            min_metallicity_by_thousand=int(min_metallicity * 1e3),
            max_metallicity_by_thousand=int(max_metallicity * 1e3)))

    spline = linear_estimation(x=(min_metallicity, max_metallicity),
                               y=(min_luminosity, max_luminosity))
    luminosity = spline(star_metallicity)

    spline = linear_estimation(x=(min_metallicity, max_metallicity),
                               y=(min_effective_temperature,
                                  max_effective_temperature))
    effective_temperature = spline(star_metallicity)

    min_mass_index = find_mass_index(star_mass=star_mass,
                                     mass_grid=mass_grid)
    rows_count_1 = rows_counts[min_mass_index]
    rows_count_2 = rows_counts[min_mass_index + 1]

    estimated_magnitude = partial(get_color_magnitude,
                                  star_mass=star_mass,
                                  star_luminosity=star_luminosity,
                                  rows_count_1=rows_count_1,
                                  rows_count_2=rows_count_2,
                                  min_mass_index=min_mass_index,
                                  luminosity_grid=luminosity_grid,
                                  mass_grid=mass_grid)

    u_ubvri_absolute = estimated_magnitude(
            magnitude_grid=color_table['u_ubvri_absolute'])
    b_ubvri_absolute = estimated_magnitude(
            magnitude_grid=color_table['b_ubvri_absolute'])
    v_ubvri_absolute = estimated_magnitude(
            magnitude_grid=color_table['v_ubvri_absolute'])
    r_ubvri_absolute = estimated_magnitude(
            magnitude_grid=color_table['r_ubvri_absolute'])
    i_ubvri_absolute = estimated_magnitude(
            magnitude_grid=color_table['i_ubvri_absolute'])

    star['luminosity'] = -luminosity
    star['effective_temperature'] = effective_temperature
    star['u_ubvri_absolute'] = u_ubvri_absolute
    star['b_ubvri_absolute'] = b_ubvri_absolute
    star['v_ubvri_absolute'] = v_ubvri_absolute
    star['r_ubvri_absolute'] = r_ubvri_absolute
    star['i_ubvri_absolute'] = i_ubvri_absolute


def get_luminosity_effective_temperature_limits(
        *,
        star_mass: float,
        star_cooling_time: float,
        cooling_sequences: Dict[int, Dict[str, np.ndarray]],
        min_metallicity_by_thousand: int,
        max_metallicity_by_thousand: int) -> Tuple[float, ...]:
    min_metallicity_sequences = cooling_sequences[min_metallicity_by_thousand]
    max_metallicity_sequences = cooling_sequences[max_metallicity_by_thousand]

    min_luminosity = interpolate(
        star_mass=star_mass,
        star_cooling_time=star_cooling_time,
        mass_grid=min_metallicity_sequences['mass'],
        cooling_time_grid=min_metallicity_sequences['cooling_time'],
        pre_wd_lifetime_grid=min_metallicity_sequences['pre_wd_lifetime'],
        interest_sequence_grid=min_metallicity_sequences['luminosity'],
        rows_counts=min_metallicity_sequences['rows_counts'],
        by_logarithm=False)
    min_effective_temperature = interpolate(
        star_mass=star_mass,
        star_cooling_time=star_cooling_time,
        mass_grid=min_metallicity_sequences['mass'],
        cooling_time_grid=min_metallicity_sequences['cooling_time'],
        pre_wd_lifetime_grid=min_metallicity_sequences['pre_wd_lifetime'],
        interest_sequence_grid=min_metallicity_sequences[
            'effective_temperature'],
        rows_counts=min_metallicity_sequences['rows_counts'],
        by_logarithm=True)
    max_luminosity = interpolate(
        star_mass=star_mass,
        star_cooling_time=star_cooling_time,
        mass_grid=max_metallicity_sequences['mass'],
        cooling_time_grid=max_metallicity_sequences['cooling_time'],
        pre_wd_lifetime_grid=max_metallicity_sequences['pre_wd_lifetime'],
        interest_sequence_grid=max_metallicity_sequences['luminosity'],
        rows_counts=max_metallicity_sequences['rows_counts'],
        by_logarithm=False)
    max_effective_temperature = interpolate(
        star_mass=star_mass,
        star_cooling_time=star_cooling_time,
        mass_grid=max_metallicity_sequences['mass'],
        cooling_time_grid=max_metallicity_sequences['cooling_time'],
        pre_wd_lifetime_grid=max_metallicity_sequences['pre_wd_lifetime'],
        interest_sequence_grid=max_metallicity_sequences[
            'effective_temperature'],
        rows_counts=max_metallicity_sequences['rows_counts'],
        by_logarithm=True)

    return (min_luminosity,
            max_luminosity,
            min_effective_temperature,
            max_effective_temperature)


def interpolate(*,
                star_mass: float,
                star_cooling_time: float,
                mass_grid: np.ndarray,
                cooling_time_grid: np.ndarray,
                pre_wd_lifetime_grid: np.ndarray,
                interest_sequence_grid: np.ndarray,
                rows_counts: np.ndarray,
                by_logarithm: bool,
                one_model: bool = False) -> float:
    extrapolate_by_mass_partial = partial(
            extrapolate_by_mass,
            star_cooling_time=star_cooling_time,
            star_mass=star_mass,
            mass_grid=mass_grid,
            cooling_time_grid=cooling_time_grid,
            pre_wd_lifetime_grid=pre_wd_lifetime_grid,
            interest_sequence_grid=interest_sequence_grid,
            rows_counts=rows_counts,
            by_logarithm=by_logarithm,
            one_model=one_model)

    if star_mass < mass_grid[0]:
        return extrapolate_by_mass_partial(min_mass_index=0)

    if star_mass >= mass_grid[-1]:
        return extrapolate_by_mass_partial(
                min_mass_index=mass_grid.size() - 1)

    return interpolate_by_mass(
        star_mass=star_mass,
        star_cooling_time=star_cooling_time,
        mass_grid=mass_grid,
        cooling_time_grid=cooling_time_grid,
        pre_wd_lifetime_grid=pre_wd_lifetime_grid,
        interest_sequence_grid=interest_sequence_grid,
        rows_counts=rows_counts,
        by_logarithm=by_logarithm,
        one_model=one_model)


def extrapolate_by_mass(*,
                        star_cooling_time: float,
                        star_mass: float,
                        min_mass_index: int,
                        mass_grid: np.ndarray,
                        cooling_time_grid: np.ndarray,
                        pre_wd_lifetime_grid: np.ndarray,
                        interest_sequence_grid: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    xm = partial(get_xm,
                 star_cooling_time=star_cooling_time,
                 cooling_time_grid=cooling_time_grid,
                 pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                 interest_sequence_grid=interest_sequence_grid,
                 rows_counts=rows_counts,
                 by_logarithm=by_logarithm,
                 one_model=one_model)
    xm1 = xm(mass_index=min_mass_index)
    xm2 = xm(mass_index=min_mass_index + 1)

    min_mass = mass_grid[min_mass_index]
    max_mass = mass_grid[min_mass_index + 1]

    spline = linear_estimation(x=(min_mass, max_mass),
                               y=(xm1, xm2))
    return spline(star_mass)


def get_mass_index(*,
                   star_mass: float,
                   mass_grid: np.ndarray) -> int:
    for row_index in range(mass_grid.size() - 1):
        if mass_grid[row_index] <= star_mass < mass_grid[row_index + 1]:
            return row_index


def interpolate_by_mass(*,
                        star_mass: float,
                        star_cooling_time: float,
                        mass_grid: np.ndarray,
                        cooling_time_grid: np.ndarray,
                        pre_wd_lifetime_grid: np.ndarray,
                        interest_sequence_grid: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    min_mass_index = get_mass_index(star_mass=star_mass,
                                    mass_grid=mass_grid)
    max_mass_index = min_mass_index + 1

    extrapolated_xm = partial(get_extrapolated_xm,
                              star_cooling_time=star_cooling_time,
                              cooling_time_grid=cooling_time_grid,
                              pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                              interest_sequence_grid=interest_sequence_grid,
                              by_logarithm=by_logarithm,
                              one_model=one_model)

    if star_cooling_time < cooling_time_grid[min_mass_index, 0]:
        x1 = extrapolated_xm(min_row_index=1,
                             mass_index=min_mass_index)
        case_1 = 1
    elif star_cooling_time >= cooling_time_grid[min_mass_index,
                                                rows_counts[min_mass_index]]:
        rows_count = rows_counts[min_mass_index]
        x1 = extrapolated_xm(min_row_index=rows_count,
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
        x3 = extrapolated_xm(min_row_index=1,
                             mass_index=max_mass_index)
        case_2 = 1
    elif star_cooling_time >= cooling_time_grid[max_mass_index,
                                                rows_counts[max_mass_index]]:
        rows_count = rows_counts[max_mass_index]
        x3 = extrapolated_xm(min_row_index=rows_count,
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

    min_mass = mass_grid[min_mass_index]
    max_mass = mass_grid[max_mass_index]
    get_linear_estimation = partial(linear_estimation,
                                    x=(min_mass, max_mass))

    if case_1 == 0 and case_2 == 0:
        spline = get_linear_estimation(y=(y1, y3))
        ym1 = spline(star_mass)

        spline = get_linear_estimation(y=(y2, y4))
        ym2 = spline(star_mass)

        spline = get_linear_estimation(y=(x1, x3))
        xm1 = spline(star_mass)

        spline = get_linear_estimation(y=(x2, x4))
        xm2 = spline(star_mass)

        spline = linear_estimation(x=(ym1, ym2),
                                   y=(xm1, xm2))
        return spline(star_cooling_time)

    if case_1 == 0 and case_2 == 1:
        spline = linear_estimation(x=(y1, y2),
                                   y=(x1, x2))
        xm1 = spline(star_cooling_time)

        spline = get_linear_estimation(y=(xm1, x3))
        return spline(star_mass)

    if case_1 == 1 and case_2 == 0:
        spline = linear_estimation(x=(y3, y4),
                                   y=(x3, x4))
        xm2 = spline(star_cooling_time)

        spline = get_linear_estimation(y=(x1, xm2))
        return spline(star_mass)

    spline = get_linear_estimation(y=(x1, x3))
    return spline(star_mass)


def get_xm(*,
           star_cooling_time: float,
           cooling_time_grid: np.ndarray,
           pre_wd_lifetime_grid: np.ndarray,
           interest_sequence_grid: np.ndarray,
           rows_counts: np.ndarray,
           mass_index: int,
           by_logarithm: bool,
           one_model: bool = False) -> float:
    extrapolated_xm = partial(get_extrapolated_xm,
                              star_cooling_time=star_cooling_time,
                              cooling_time_grid=cooling_time_grid,
                              pre_wd_lifetime_grid=pre_wd_lifetime_grid,
                              interest_sequence_grid=interest_sequence_grid,
                              mass_index=mass_index,
                              by_logarithm=by_logarithm,
                              one_model=one_model)

    if star_cooling_time < cooling_time_grid[mass_index, 0]:
        return extrapolated_xm(min_row_index=0)

    rows_count = rows_counts[mass_index]

    if star_cooling_time > cooling_time_grid[mass_index, rows_count]:
        return extrapolated_xm(min_row_index=rows_count - 1)

    for row_index in range(rows_count - 1):
        if (cooling_time_grid[mass_index, row_index] <= star_cooling_time
                <= cooling_time_grid[mass_index, row_index + 1]):
            return estimate_at(
                    star_cooling_time,
                    x=(cooling_time_grid[mass_index, row_index],
                       cooling_time_grid[mass_index, row_index + 1]),
                    y=(interest_sequence_grid[mass_index, row_index],
                       interest_sequence_grid[mass_index, row_index + 1]))


def get_extrapolated_xm(*,
                        star_cooling_time: float,
                        cooling_time_grid: np.ndarray,
                        pre_wd_lifetime_grid: np.ndarray,
                        interest_sequence_grid: np.ndarray,
                        mass_index: int,
                        min_row_index: int,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    if one_model:
        return estimate_at(
                star_cooling_time,
                x=(cooling_time_grid[mass_index, min_row_index],
                   cooling_time_grid[mass_index, min_row_index + 1]),
                y=(interest_sequence_grid[mass_index, min_row_index],
                   interest_sequence_grid[mass_index, min_row_index + 1]))
    if by_logarithm:
        return 10. ** estimate_at(
                star_cooling_time,
                x=(log10(cooling_time_grid[mass_index, min_row_index]
                         + pre_wd_lifetime_grid[mass_index]),
                   log10(cooling_time_grid[mass_index, min_row_index + 1]
                         + pre_wd_lifetime_grid[mass_index])),
                y=(log10(interest_sequence_grid[mass_index,
                                                min_row_index]),
                   log10(interest_sequence_grid[mass_index,
                                                min_row_index + 1])))

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
    for mass_index in range(mass_grid.size() - 1):
        if (mass_grid[mass_index]
                < star_mass
                <= mass_grid[mass_index + 1]):
            return mass_index


def get_color_magnitude(*,
                        star_mass: float,
                        star_luminosity: float,
                        rows_count_1: int,
                        rows_count_2: int,
                        min_mass_index: int,
                        luminosity_grid: np.ndarray,
                        magnitude_grid: np.ndarray,
                        mass_grid: np.ndarray,
                        extrapolating_case: bool = False) -> float:
    if (star_luminosity > luminosity_grid[min_mass_index, 0]
            or star_luminosity > luminosity_grid[min_mass_index + 1, 0]):
        extrapolating_case = True
        row_index_1 = 0
        row_index_2 = 0

    if (star_luminosity < luminosity_grid[min_mass_index, rows_count_1]
            or star_luminosity < luminosity_grid[min_mass_index + 1,
                                                 rows_count_2]):
        extrapolating_case = True
        row_index_1 = rows_count_1
        row_index_2 = rows_count_2

    if extrapolating_case:
        min_magnitude = estimate_at(
                star_luminosity,
                x=(luminosity_grid[min_mass_index, row_index_1],
                   luminosity_grid[min_mass_index, row_index_1 + 1]),
                y=(magnitude_grid[min_mass_index, row_index_1],
                   magnitude_grid[min_mass_index, row_index_1 + 1]))
        max_magnitude = estimate_at(
                star_luminosity,
                x=(luminosity_grid[min_mass_index + 1, row_index_2],
                   luminosity_grid[min_mass_index + 1, row_index_2 + 1]),
                y=(magnitude_grid[min_mass_index + 1, row_index_2],
                   magnitude_grid[min_mass_index + 1, row_index_2 + 1]))
        magnitude = estimate_at(star_mass,
                                x=(mass_grid[min_mass_index],
                                   mass_grid[min_mass_index + 1]),
                                y=(min_magnitude, max_magnitude))
        return max(0, magnitude)

    find_row_index = partial(find_index,
                             luminosity=star_luminosity,
                             luminosity_grid=luminosity_grid)
    row_index_1 = find_row_index(rows_count=rows_count_1,
                                 mass_index=min_mass_index)
    row_index_2 = find_row_index(rows_count=rows_count_2,
                                 mass_index=min_mass_index + 1)
    min_magnitude = estimate_at(
            star_luminosity,
            x=(luminosity_grid[min_mass_index, row_index_1],
               luminosity_grid[min_mass_index, row_index_1 + 1]),
            y=(magnitude_grid[min_mass_index, row_index_1],
               magnitude_grid[min_mass_index, row_index_1 + 1]))
    max_magnitude = estimate_at(
            star_luminosity,
            x=(luminosity_grid[min_mass_index + 1, row_index_2],
               luminosity_grid[min_mass_index + 1, row_index_2 + 1]),
            y=(magnitude_grid[min_mass_index + 1, row_index_2],
               magnitude_grid[min_mass_index + 1, row_index_2 + 1]))

    return estimate_at(star_mass,
                       x=(mass_grid[0], mass_grid[1]),
                       y=(min_magnitude, max_magnitude))


def find_index(*,
               rows_count: int,
               star_luminosity: float,
               luminosity_grid: np.ndarray,
               mass_index: int) -> int:
    for row_index in range(rows_count - 1):
        if (luminosity_grid[mass_index, row_index + 1]
                <= star_luminosity
                <= luminosity_grid[mass_index, row_index]):
            return row_index


def estimate_at(x0: float,
                *,
                x: Tuple[float, float],
                y: Tuple[float, float]) -> float:
    spline = linear_estimation(x=x,
                               y=y)
    return spline(x0)

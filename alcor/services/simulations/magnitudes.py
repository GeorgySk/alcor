import enum
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
    carbon_oxygen_white_dwarfs_mask = (stars['mass']
                                       < max_carbon_oxygen_core_wd_mass)
    carbon_oxygen_white_dwarfs = stars[carbon_oxygen_white_dwarfs_mask]
    oxygen_neon_white_dwarfs = stars[~carbon_oxygen_white_dwarfs_mask]

    for _, star in carbon_oxygen_white_dwarfs.iterrows():
        if get_spectral_type(db_to_da_fraction) == SpectralType.DA:
            star['spectral_type'] = SpectralType.DA
            (luminosity,
             effective_temperature,
             u_ubvri_absolute,
             b_ubvri_absolute,
             v_ubvri_absolute,
             r_ubvri_absolute,
             i_ubvri_absolute) = da_db_interpolation(
                star=star,
                cooling_sequences=da_cooling_sequences,
                color_table=da_color_table,
                # TODO. can they be taken from cool.seq. keys?
                metallicities=[0.001, 0.01, 0.03, 0.06])
        else:
            star['spectral_type'] = SpectralType.DB
            (luminosity,
             effective_temperature,
             u_ubvri_absolute,
             b_ubvri_absolute,
             v_ubvri_absolute,
             r_ubvri_absolute,
             i_ubvri_absolute) = da_db_interpolation(
                star=star,
                cooling_sequences=db_cooling_sequences,
                color_table=db_color_table,
                metallicities=[0.001, 0.01, 0.06])

        star['luminosity'] = -luminosity
        star['effective_temperature'] = effective_temperature
        star['u_ubvri_absolute'] = u_ubvri_absolute
        star['b_ubvri_absolute'] = b_ubvri_absolute
        star['v_ubvri_absolute'] = v_ubvri_absolute
        star['r_ubvri_absolute'] = r_ubvri_absolute
        star['i_ubvri_absolute'] = i_ubvri_absolute

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


def get_spectral_type(db_to_da_fraction: float) -> enum.Enum:
    if np.random.rand() < db_to_da_fraction:
        return SpectralType.DB
    return SpectralType.DA


def one_interpolation(star: pd.Series,
                      color_table: Dict[str, np.ndarray],
                      one_model: bool = True,
                      by_logarithm: bool = False) -> Tuple[float, ...]:
    star['cooling_time'] = log10(star['cooling_time']) + 9.

    luminosity = interpolate(star=star,
                             cooling_or_color_sequence=color_table,
                             interest_sequence='luminosity',
                             by_logarithm=by_logarithm,
                             one_model=one_model)
    v_ubvri_absolute = interpolate(star=star,
                                   cooling_or_color_sequence=color_table,
                                   interest_sequence='v_ubvri_absolute',
                                   by_logarithm=by_logarithm,
                                   one_model=one_model)
    bv_ubvri = interpolate(star=star,
                           cooling_or_color_sequence=color_table,
                           interest_sequence='bv_ubvri',
                           by_logarithm=by_logarithm,
                           one_model=one_model)
    vi_ubvri = interpolate(star=star,
                           cooling_or_color_sequence=color_table,
                           interest_sequence='vi_ubvri',
                           by_logarithm=by_logarithm,
                           one_model=one_model)
    vr_ubvri = interpolate(star=star,
                           cooling_or_color_sequence=color_table,
                           interest_sequence='vr_ubvri',
                           by_logarithm=by_logarithm,
                           one_model=one_model)
    uv_ubvri = interpolate(star=star,
                           cooling_or_color_sequence=color_table,
                           interest_sequence='uv_ubvri',
                           by_logarithm=by_logarithm,
                           one_model=one_model)
    log_effective_temperature = interpolate(
        star=star,
        cooling_or_color_sequence=color_table,
        interest_sequence='log_effective_temperature',
        by_logarithm=by_logarithm,
        one_model=one_model)

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


def da_db_interpolation(star: pd.Series,
                        cooling_sequences: Dict[int, Dict[str, np.ndarray]],
                        color_table: Dict[str, np.ndarray],
                        metallicities: List[float]) -> Tuple[float, ...]:
    for metallicity_index in range(len(metallicities) - 1):
        if (metallicities[metallicity_index] <= star['metallicity']
                < metallicities[metallicity_index + 1]):
            min_metallicity = metallicities[metallicity_index]
            max_metallicity = metallicities[metallicity_index + 1]

            (min_luminosity,
             max_luminosity,
             min_effective_temperature,
             max_effective_temperature) = (
                get_luminosity_effective_temperature_limits(
                    star=star,
                    cooling_sequences=cooling_sequences,
                    min_metallicity_by_thousand=int(min_metallicity * 1e3),
                    max_metallicity_by_thousand=int(max_metallicity * 1e3)))
            break

    # TODO: this looks like linear extrapolation. implement function
    luminosity = (min_luminosity
                  + (max_luminosity - min_luminosity)
                    * (star['metallicity'] - min_metallicity)
                    / (max_metallicity - min_metallicity))
    effective_temperature = (min_effective_temperature
                             + (max_effective_temperature
                                - min_effective_temperature)
                               * (star['metallicity'] - min_metallicity)
                               / (max_metallicity - min_metallicity))

    (u_ubvri_absolute,
     b_ubvri_absolute,
     v_ubvri_absolute,
     r_ubvri_absolute,
     i_ubvri_absolute) = interpolate_magnitudes(star_mass=star,
                                                color_table=color_table,
                                                luminosity=luminosity)

    return (luminosity,
            effective_temperature,
            u_ubvri_absolute,
            b_ubvri_absolute,
            v_ubvri_absolute,
            r_ubvri_absolute,
            i_ubvri_absolute)


def get_luminosity_effective_temperature_limits(
        *,
        star: pd.Series,
        cooling_sequences: Dict[int, Dict[str, np.ndarray]],
        min_metallicity_by_thousand: int,
        max_metallicity_by_thousand: int) -> Tuple[float, ...]:
    min_metallicity_sequences = cooling_sequences[min_metallicity_by_thousand]
    max_metallicity_sequences = cooling_sequences[max_metallicity_by_thousand]

    min_luminosity = interpolate(
        star=star,
        cooling_or_color_sequence=min_metallicity_sequences,
        interest_sequence='luminosity',
        by_logarithm=False)
    min_effective_temperature = interpolate(
        star=star,
        cooling_or_color_sequence=min_metallicity_sequences,
        interest_sequence='effective_temperature',
        by_logarithm=True)
    max_luminosity = interpolate(
        star=star,
        cooling_or_color_sequence=max_metallicity_sequences,
        interest_sequence='luminosity',
        by_logarithm=False)
    max_effective_temperature = interpolate(
        star=star,
        cooling_or_color_sequence=max_metallicity_sequences,
        interest_sequence='effective_temperature',
        by_logarithm=True)

    return (min_luminosity,
            max_luminosity,
            min_effective_temperature,
            max_effective_temperature)


def interpolate(*,
                star: pd.Series,
                cooling_or_color_sequence: Dict[str, np.ndarray],
                interest_sequence: str,
                by_logarithm: bool,
                one_model: bool = False) -> float:
    grid_masses = cooling_or_color_sequence['mass']
    grid_cooling_times = cooling_or_color_sequence['cooling_time']
    grid_pre_wd_lifetimes = cooling_or_color_sequence['pre_wd_lifetime']
    rows_counts = cooling_or_color_sequence['rows_counts']
    grid_interest_sequence = cooling_or_color_sequence[interest_sequence]

    extrapolate_by_mass_partial = partial(
            extrapolate_by_mass,
            star=star,
            mass=grid_masses,
            cooling_time=grid_cooling_times,
            pre_wd_lifetime=grid_pre_wd_lifetimes,
            interest_sequence=grid_interest_sequence,
            rows_counts=rows_counts,
            by_logarithm=by_logarithm,
            one_model=one_model)

    star_mass = star['mass']

    if star_mass < grid_masses[0]:
        return extrapolate_by_mass_partial(min_mass_index=0)

    if star_mass >= grid_masses[-1]:
        return extrapolate_by_mass_partial(
                min_mass_index=grid_masses.size() - 1)

    return interpolate_by_mass(
        star=star,
        mass=grid_masses,
        cooling_time=grid_cooling_times,
        pre_wd_lifetime=grid_pre_wd_lifetimes,
        interest_sequence=grid_interest_sequence,
        rows_counts=rows_counts,
        by_logarithm=by_logarithm,
        one_model=one_model)


def extrapolate_by_mass(*,
                        star: pd.Series,
                        min_mass_index: int,
                        mass: np.ndarray,
                        cooling_time: np.ndarray,
                        pre_wd_lifetime: np.ndarray,
                        interest_sequence: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    xm = partial(get_xm,
                 star_cooling_time=star['cooling_time'],
                 cooling_time=cooling_time,
                 pre_wd_lifetime=pre_wd_lifetime,
                 interest_sequence=interest_sequence,
                 rows_counts=rows_counts,
                 by_logarithm=by_logarithm,
                 one_model=one_model)
    xm1 = xm(mass_index=min_mass_index)
    xm2 = xm(mass_index=min_mass_index + 1)

    min_mass = mass[min_mass_index]
    max_mass = mass[min_mass_index + 1]

    spline = linear_estimation(x=(min_mass, max_mass),
                               y=(xm1, xm2))
    return spline(star['mass'])


def get_mass_index(*,
                   star_mass: float,
                   grid_masses: np.ndarray) -> int:
    for row_index in range(grid_masses.size() - 1):
        if grid_masses[row_index] <= star_mass < grid_masses[row_index + 1]:
            return row_index


def interpolate_by_mass(*,
                        star: pd.Series,
                        mass: np.ndarray,
                        cooling_time: np.ndarray,
                        pre_wd_lifetime: np.ndarray,
                        interest_sequence: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    star_mass = star['mass']
    star_cooling_time = star['cooling_time']

    min_mass_index = get_mass_index(star_mass=star_mass,
                                    grid_masses=mass)
    max_mass_index = min_mass_index + 1

    extrapolated_xm = partial(get_extrapolated_xm,
                              star_cooling_time=star_cooling_time,
                              cooling_time=cooling_time,
                              pre_wd_lifetime=pre_wd_lifetime,
                              interest_sequence=interest_sequence,
                              by_logarithm=by_logarithm,
                              one_model=one_model)

    if star_cooling_time < cooling_time[min_mass_index, 0]:
        x1 = extrapolated_xm(min_row_index=1,
                             mass_index=min_mass_index)
        case_1 = 1
    elif star_cooling_time >= cooling_time[min_mass_index,
                                           rows_counts[min_mass_index]]:
        rows_count = rows_counts[min_mass_index]
        x1 = extrapolated_xm(min_row_index=rows_count,
                             mass_index=min_mass_index)
        case_1 = 1
    else:
        for row_index in range(rows_counts[min_mass_index] - 1):
            if (cooling_time[min_mass_index, row_index]
                    <= star_cooling_time
                    <= cooling_time[min_mass_index, row_index + 1]):
                y1 = cooling_time[min_mass_index, row_index]
                y2 = cooling_time[min_mass_index, row_index + 1]
                x1 = interest_sequence[min_mass_index, row_index]
                x2 = interest_sequence[min_mass_index, row_index + 1]
                case_1 = 0

    if star_cooling_time < cooling_time[max_mass_index, 0]:
        x3 = extrapolated_xm(min_row_index=1,
                             mass_index=max_mass_index)
        case_2 = 1
    elif star_cooling_time >= cooling_time[max_mass_index,
                                           rows_counts[max_mass_index]]:
        rows_count = rows_counts[max_mass_index]
        x3 = extrapolated_xm(min_row_index=rows_count,
                             mass_index=max_mass_index)
        case_2 = 1
    else:
        for row_index in range(rows_counts[max_mass_index] - 1):
            if (cooling_time[max_mass_index, row_index]
                    <= star_cooling_time
                    <= cooling_time[max_mass_index, row_index + 1]):
                y3 = cooling_time[max_mass_index, row_index]
                y4 = cooling_time[max_mass_index, row_index + 1]
                x3 = interest_sequence[max_mass_index, row_index]
                x4 = interest_sequence[max_mass_index, row_index + 1]
                case_2 = 0

    min_mass = mass[min_mass_index]
    max_mass = mass[max_mass_index]
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
           cooling_time: np.ndarray,
           pre_wd_lifetime: np.ndarray,
           interest_sequence: np.ndarray,
           rows_counts: np.ndarray,
           mass_index: int,
           by_logarithm: bool,
           one_model: bool = False) -> float:
    extrapolated_xm = partial(get_extrapolated_xm,
                              star_cooling_time=star_cooling_time,
                              cooling_time=cooling_time,
                              pre_wd_lifetime=pre_wd_lifetime,
                              interest_sequence=interest_sequence,
                              mass_index=mass_index,
                              by_logarithm=by_logarithm,
                              one_model=one_model)

    if star_cooling_time < cooling_time[mass_index, 0]:
        return extrapolated_xm(min_row_index=0)

    rows_count = rows_counts[mass_index]

    if star_cooling_time > cooling_time[mass_index, rows_count]:
        return extrapolated_xm(min_row_index=rows_count - 1)

    for row_index in range(rows_count - 1):
        if (cooling_time[mass_index, row_index] <= star_cooling_time
                <= cooling_time[mass_index, row_index + 1]):
            spline = linear_estimation(
                    x=(cooling_time[mass_index, row_index],
                       cooling_time[mass_index, row_index + 1]),
                    y=(interest_sequence[mass_index, row_index],
                       interest_sequence[mass_index, row_index + 1]))
            return spline(star_cooling_time)


def get_extrapolated_xm(*,
                        star_cooling_time: float,
                        cooling_time: np.ndarray,
                        pre_wd_lifetime: np.ndarray,
                        interest_sequence: np.ndarray,
                        mass_index: int,
                        min_row_index: int,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    if one_model:
        spline = linear_estimation(
                x=(cooling_time[mass_index, min_row_index],
                   cooling_time[mass_index, min_row_index + 1]),
                y=(interest_sequence[mass_index, min_row_index],
                   interest_sequence[mass_index, min_row_index + 1]))
        return spline(star_cooling_time)

    x0 = log10(cooling_time[mass_index, min_row_index]
               + pre_wd_lifetime[mass_index])
    x1 = log10(cooling_time[mass_index, min_row_index + 1]
               + pre_wd_lifetime[mass_index])

    if by_logarithm:
        y0 = log10(interest_sequence[mass_index, min_row_index])
        y1 = log10(interest_sequence[mass_index, min_row_index + 1])
        spline = linear_estimation(x=(x0, x1),
                                   y=(y0, y1))
        return 10.0 ** spline(star_cooling_time)

    spline = linear_estimation(
            x=(x0, x1),
            y=(interest_sequence[mass_index, min_row_index],
               interest_sequence[mass_index, min_row_index + 1]))
    return spline(star_cooling_time)


def interpolate_magnitudes(*,
                           star_mass: float,
                           color_table: Dict[str, np.ndarray],
                           luminosity: float) -> Tuple[float, ...]:
    grid_masses = color_table['mass']
    rows_counts = color_table['rows_counts']

    if star_mass <= grid_masses[0]:
        min_mass_index = 0
    elif star_mass > grid_masses[-1]:
        # Index of element before the last one
        min_mass_index = -2
    else:
        min_mass_index = find_mass_index(star_mass=star_mass,
                                         grid_masses=grid_masses)

    rows_count_1 = rows_counts[min_mass_index]
    rows_count_2 = rows_counts[min_mass_index + 1]

    return get_magnitudes(star_mass=star_mass,
                          luminosity=luminosity,
                          color_table=color_table,
                          rows_count_1=rows_count_1,
                          rows_count_2=rows_count_2,
                          min_mass_index=min_mass_index)


def find_mass_index(*,
                    star_mass: float,
                    grid_masses: np.ndarray) -> int:
    for mass_index in range(grid_masses.size() - 1):
        if (grid_masses[mass_index]
                < star_mass
                <= grid_masses[mass_index + 1]):
            return mass_index


def get_magnitudes(*,
                   star_mass: float,
                   luminosity: float,
                   color_table: Dict[str, np.ndarray],
                   rows_count_1: int,
                   rows_count_2: int,
                   min_mass_index: int) -> Tuple[float, ...]:
    luminosity_grid = color_table['luminosity']

    extrapolated_magnitudes = partial(
            get_extrapolated_magnitudes_by_luminosity,
            star_mass=star_mass,
            luminosity=luminosity,
            color_table=color_table,
            mass_index=min_mass_index)

    if (luminosity > luminosity_grid[min_mass_index, 0]
            or luminosity > luminosity_grid[min_mass_index + 1, 0]):
        return extrapolated_magnitudes(row_index_1=0,
                                       row_index_2=0)

    if (luminosity < luminosity_grid[min_mass_index, rows_count_1]
            or luminosity < luminosity_grid[min_mass_index + 1, rows_count_2]):
        return extrapolated_magnitudes(row_index_1=rows_count_1,
                                       row_index_2=rows_count_2)

    return get_interpolated_magnitudes_by_luminosity(
        star_mass=star_mass,
        rows_count_1=rows_count_1,
        rows_count_2=rows_count_2,
        color_table=color_table,
        luminosity=luminosity,
        mass_index=min_mass_index)


def get_interpolated_magnitudes_by_luminosity(
        *,
        star_mass: float,
        rows_count_1: int,
        rows_count_2: int,
        color_table: Dict[str, np.ndarray],
        luminosity: float,
        mass_index: int) -> Tuple[float, ...]:
    luminosity_grid = color_table['luminosity']

    find_row_index = partial(find_index,
                             luminosity=luminosity,
                             luminosity_grid=luminosity_grid)
    row_index_1 = find_row_index(rows_count=rows_count_1,
                                 mass_index=mass_index)
    row_index_2 = find_row_index(rows_count=rows_count_2,
                                 mass_index=mass_index + 1)

    get_magnitude = partial(get_interpolated_magnitude,
                            star_mass=star_mass,
                            star_luminosity=luminosity,
                            table_luminosity=luminosity_grid,
                            table_mass=color_table['mass'],
                            row_index_1=row_index_1,
                            row_index_2=row_index_2,
                            mass_index=mass_index)

    return (get_magnitude(table_magnitude=color_table['u_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['b_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['v_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['r_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['i_ubvri_absolute']))


def find_index(*,
               rows_count: int,
               luminosity: float,
               luminosity_grid: np.ndarray,
               mass_index: int) -> int:
    for row_index in range(rows_count - 1):
        if (luminosity_grid[mass_index, row_index + 1]
                <= luminosity
                <= luminosity_grid[mass_index, row_index]):
            return row_index


def get_extrapolated_magnitudes_by_luminosity(
        *,
        star_mass: float,
        luminosity: float,
        color_table: Dict[str, np.ndarray],
        row_index_1: int,
        row_index_2: int,
        mass_index: int) -> Tuple[float, ...]:
    get_magnitude = partial(get_abs_magnitude,
                            star_mass=star_mass,
                            star_luminosity=luminosity,
                            table_luminosity=color_table['luminosity'],
                            table_mass=color_table['mass'],
                            row_index_1=row_index_1,
                            row_index_2=row_index_2,
                            mass_index=mass_index)

    return (get_magnitude(table_magnitude=color_table['u_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['b_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['v_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['r_ubvri_absolute']),
            get_magnitude(table_magnitude=color_table['i_ubvri_absolute']))


# TODO: looks more like extrapolated
def get_interpolated_magnitude(*,
                               star_mass: float,
                               star_luminosity: float,
                               table_luminosity: np.ndarray,
                               table_magnitude: np.ndarray,
                               table_mass: np.ndarray,
                               row_index_1: int,
                               row_index_2: int,
                               mass_index: int) -> float:
    c_1_spline = linear_estimation(
            x=(table_luminosity[mass_index, row_index_1],
               table_luminosity[mass_index, row_index_1 + 1]),
            y=(table_magnitude[mass_index, row_index_1],
               table_magnitude[mass_index, row_index_1 + 1]))
    c_1 = c_1_spline(star_luminosity)

    c_2_spline = linear_estimation(
            x=(table_luminosity[mass_index + 1, row_index_2],
               table_luminosity[mass_index + 1, row_index_2 + 1]),
            y=(table_magnitude[mass_index + 1, row_index_2],
               table_magnitude[mass_index + 1, row_index_2 + 1]))
    c_2 = c_2_spline(star_luminosity)

    magnitude_spline = linear_estimation(x=(table_mass[0], table_mass[1]),
                                         y=(c_1, c_2))
    return magnitude_spline(star_mass)


# TODO: is this interpolation?
def get_abs_magnitude(*,
                      star_mass: float,
                      star_luminosity: float,
                      table_magnitude: np.ndarray,
                      table_luminosity: np.ndarray,
                      table_mass: np.ndarray,
                      row_index_1: int,
                      row_index_2: int,
                      mass_index: int) -> float:
    c_1_spline = linear_estimation(
            x=(table_luminosity[mass_index, row_index_1],
               table_luminosity[mass_index, row_index_1 + 1]),
            y=(table_magnitude[mass_index, row_index_1],
               table_magnitude[mass_index, row_index_1 + 1]))
    c_1 = c_1_spline(star_luminosity)

    c_2_spline = linear_estimation(
            x=(table_luminosity[mass_index + 1, row_index_2],
               table_luminosity[mass_index + 1, row_index_2 + 1]),
            y=(table_magnitude[mass_index + 1, row_index_2],
               table_magnitude[mass_index + 1, row_index_2 + 1]))
    c_2 = c_2_spline(star_luminosity)

    magnitude_spline = linear_estimation(
            x=(table_mass[mass_index], table_mass[mass_index + 1]),
            y=(c_1, c_2))
    abs_magnitude = magnitude_spline(star_mass)

    return max(0, abs_magnitude)

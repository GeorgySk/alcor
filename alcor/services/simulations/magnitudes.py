import enum
from math import log10
from typing import (Dict,
                    Tuple,
                    List)

import numpy as np

from alcor.models import Star
from alcor.models.star import SpectralTypeEnum
from alcor.types import CoolingSequencesType

GRAVITATIONAL_CONST_CM_S_KG = 6.67e-5
SOLAR_MASS_KG = 1.989e30


def assign_magnitudes(*,
                      stars: List[Star],
                      max_carbon_oxygen_core_wd_mass: float = 1.14,
                      db_to_da_fraction: float = 0.2,
                      da_cooling_sequences: CoolingSequencesType,
                      da_color_table: Dict[str, np.ndarray],
                      db_cooling_sequences: CoolingSequencesType,
                      db_color_table: Dict[str, np.ndarray],
                      one_color_table: Dict[str, np.ndarray]
                      ) -> List[Star]:
    co_wds = [star
              for star in stars
              if star.mass < max_carbon_oxygen_core_wd_mass]
    one_wds = [star
               for star in stars
               if star.mass >= max_carbon_oxygen_core_wd_mass]

    for star in co_wds:
        if get_spectral_type(db_to_da_fraction) == SpectralTypeEnum.DA:
            star.spectral_type = SpectralTypeEnum.DA
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
            star.spectral_type = SpectralTypeEnum.DB
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

        star.luminosity = -luminosity
        star.effective_temperature = effective_temperature
        star.u_ubvri_absolute = u_ubvri_absolute
        star.b_ubvri_absolute = b_ubvri_absolute
        star.v_ubvri_absolute = v_ubvri_absolute
        star.r_ubvri_absolute = r_ubvri_absolute
        star.i_ubvri_absolute = i_ubvri_absolute

    for star in one_wds:
        star.spectral_type = SpectralTypeEnum.ONe
        (luminosity,
         effective_temperature,
         u_ubvri_absolute,
         b_ubvri_absolute,
         v_ubvri_absolute,
         r_ubvri_absolute,
         i_ubvri_absolute) = one_interpolation(star=star,
                                               color_table=one_color_table)

        star.luminosity = -luminosity
        star.effective_temperature = effective_temperature
        star.u_ubvri_absolute = u_ubvri_absolute
        star.b_ubvri_absolute = b_ubvri_absolute
        star.v_ubvri_absolute = v_ubvri_absolute
        star.r_ubvri_absolute = r_ubvri_absolute
        star.i_ubvri_absolute = i_ubvri_absolute

    return co_wds + one_wds


def get_spectral_type(db_to_da_fraction: float) -> enum.Enum:
    if np.random.rand() < db_to_da_fraction:
        return SpectralTypeEnum.DB
    return SpectralTypeEnum.DA


def one_interpolation(star: Star,
                      color_table: Dict[str, np.ndarray]
                      ) -> Tuple[float, ...]:
    one_model = True
    by_logarithm = False

    star.cooling_time = log10(star.cooling_time) + 9.

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


def da_db_interpolation(star: Star,
                        cooling_sequences: CoolingSequencesType,
                        color_table: Dict[str, np.ndarray],
                        metallicities: List[float]
                        ) -> Tuple[float, ...]:
    for metallicity_index in range(len(metallicities) - 1):
        if (metallicities[metallicity_index] <= star.metallicity
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
                    min_metallicity_by_thousand=int(min_metallicity * 1000),
                    max_metallicity_by_thousand=int(max_metallicity * 1000)))
            break

    luminosity = (min_luminosity
                  + (max_luminosity - min_luminosity)
                    * (star.metallicity - min_metallicity)
                    / (max_metallicity - min_metallicity))
    effective_temperature = (min_effective_temperature
                             + (max_effective_temperature
                                - min_effective_temperature)
                               * (star.metallicity - min_metallicity)
                               / (max_metallicity - min_metallicity))

    (u_ubvri_absolute,
     b_ubvri_absolute,
     v_ubvri_absolute,
     r_ubvri_absolute,
     i_ubvri_absolute) = interpolate_magnitudes(star=star,
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
        star: Star,
        cooling_sequences: CoolingSequencesType,
        min_metallicity_by_thousand: int,
        max_metallicity_by_thousand: int) -> Tuple[float, ...]:
    min_luminosity = interpolate(
        star,
        cooling_sequences[min_metallicity_by_thousand],
        interest_sequence='luminosity',
        by_logarithm=False)
    min_effective_temperature = interpolate(
        star,
        cooling_sequences[min_metallicity_by_thousand],
        interest_sequence='effective_temperature',
        by_logarithm=True)
    max_luminosity = interpolate(
        star,
        cooling_sequences[max_metallicity_by_thousand],
        interest_sequence='luminosity',
        by_logarithm=False)
    max_effective_temperature = interpolate(
        star,
        cooling_sequences[max_metallicity_by_thousand],
        interest_sequence='effective_temperature',
        by_logarithm=True)

    return (min_luminosity, max_luminosity,
            min_effective_temperature, max_effective_temperature)


def interpolate(star: Star,
                cooling_or_color_sequence: Dict[str, np.ndarray],
                interest_sequence: str,
                by_logarithm: bool,
                one_model: bool = False) -> float:
    if star.mass < cooling_or_color_sequence['mass'][0]:
        return extrapolate_by_mass(
            star=star,
            min_mass_index=0,
            mass=cooling_or_color_sequence['mass'],
            cooling_time=cooling_or_color_sequence['cooling_time'],
            pre_wd_lifetime=cooling_or_color_sequence['pre_wd_lifetime'],
            interest_sequence=cooling_or_color_sequence[interest_sequence],
            rows_counts=cooling_or_color_sequence['rows_counts'],
            by_logarithm=by_logarithm,
            one_model=one_model)

    if star.mass >= cooling_or_color_sequence['mass'][-1]:
        return extrapolate_by_mass(
            star=star,
            min_mass_index=cooling_or_color_sequence['mass'].size() - 1,
            mass=cooling_or_color_sequence['mass'],
            cooling_time=cooling_or_color_sequence['cooling_time'],
            pre_wd_lifetime=cooling_or_color_sequence['pre_wd_lifetime'],
            interest_sequence=cooling_or_color_sequence[interest_sequence],
            rows_counts=cooling_or_color_sequence['rows_counts'],
            by_logarithm=by_logarithm,
            one_model=one_model)

    return interpolate_by_mass(
        star=star,
        mass=cooling_or_color_sequence['mass'],
        cooling_time=cooling_or_color_sequence['cooling_time'],
        pre_wd_lifetime=cooling_or_color_sequence['pre_wd_lifetime'],
        interest_sequence=cooling_or_color_sequence[interest_sequence],
        rows_counts=cooling_or_color_sequence['rows_counts'],
        by_logarithm=by_logarithm,
        one_model=one_model)


def extrapolate_by_mass(star: Star,
                        min_mass_index: int,
                        mass: np.ndarray,
                        cooling_time: np.ndarray,
                        pre_wd_lifetime: np.ndarray,
                        interest_sequence: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    min_mass = mass[min_mass_index]
    xm1 = get_xm(star=star,
                 cooling_time=cooling_time,
                 pre_wd_lifetime=pre_wd_lifetime,
                 interest_sequence=interest_sequence,
                 rows_counts=rows_counts,
                 mass_index=min_mass_index,
                 by_logarithm=by_logarithm,
                 one_model=one_model)

    max_mass = mass[min_mass_index + 1]
    xm2 = get_xm(star=star,
                 cooling_time=cooling_time,
                 pre_wd_lifetime=pre_wd_lifetime,
                 interest_sequence=interest_sequence,
                 rows_counts=rows_counts,
                 mass_index=min_mass_index + 1,
                 by_logarithm=by_logarithm,
                 one_model=one_model)

    s = (xm2 - xm1) / (max_mass - min_mass)
    t = xm2 - s * max_mass
    return s * star.mass + t


def interpolate_by_mass(star: Star,
                        mass: np.ndarray,
                        cooling_time: np.ndarray,
                        pre_wd_lifetime: np.ndarray,
                        interest_sequence: np.ndarray,
                        rows_counts: np.ndarray,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    for row_index in range(mass.size() - 1):
        if mass[row_index] <= star.mass < mass[row_index + 1]:
            min_mass_index = row_index
            max_mass_index = row_index + 1
            break

    if star.cooling_time < cooling_time[min_mass_index, 0]:
        x1 = get_extrapolated_xm(
            star=star,
            cooling_time=cooling_time,
            pre_wd_lifetime=pre_wd_lifetime,
            interest_sequence=interest_sequence,
            mass_index=min_mass_index,
            min_row_index=1,
            by_logarithm=by_logarithm,
            one_model=one_model)
        case_1 = 1
    else:
        if star.cooling_time >= cooling_time[min_mass_index,
                                             rows_counts[min_mass_index]]:
            rows_count = rows_counts[min_mass_index]
            x1 = get_extrapolated_xm(
                star=star,
                cooling_time=cooling_time,
                pre_wd_lifetime=pre_wd_lifetime,
                interest_sequence=interest_sequence,
                mass_index=min_mass_index,
                min_row_index=rows_count,
                by_logarithm=by_logarithm,
                one_model=one_model)
            case_1 = 1
        else:
            for row_index in range(rows_counts[min_mass_index] - 1):
                if (cooling_time[min_mass_index, row_index]
                        <= star.cooling_time
                        <= cooling_time[min_mass_index, row_index + 1]):
                    y1 = cooling_time[min_mass_index, row_index]
                    y2 = cooling_time[min_mass_index, row_index + 1]
                    x1 = interest_sequence[min_mass_index, row_index]
                    x2 = interest_sequence[min_mass_index, row_index + 1]
                    case_1 = 0

    if star.cooling_time < cooling_time[max_mass_index, 0]:
        x3 = get_extrapolated_xm(
            star=star,
            cooling_time=cooling_time,
            pre_wd_lifetime=pre_wd_lifetime,
            interest_sequence=interest_sequence,
            mass_index=max_mass_index,
            min_row_index=1,
            by_logarithm=by_logarithm,
            one_model=one_model)
        case_2 = 1
    else:
        if star.cooling_time >= cooling_time[max_mass_index,
                                             rows_counts[max_mass_index]]:
            rows_count = rows_counts[max_mass_index]
            x3 = get_extrapolated_xm(
                star=star,
                cooling_time=cooling_time,
                pre_wd_lifetime=pre_wd_lifetime,
                interest_sequence=interest_sequence,
                mass_index=max_mass_index,
                min_row_index=rows_count,
                by_logarithm=by_logarithm,
                one_model=one_model)
            case_2 = 1
        else:
            for row_index in range(rows_counts[max_mass_index] - 1):
                if (cooling_time[max_mass_index, row_index]
                        <= star.cooling_time
                        <= cooling_time[max_mass_index, row_index + 1]):
                    y3 = cooling_time[max_mass_index, row_index]
                    y4 = cooling_time[max_mass_index, row_index + 1]
                    x3 = interest_sequence[max_mass_index, row_index]
                    x4 = interest_sequence[max_mass_index, row_index + 1]
                    case_2 = 0

    min_mass = mass[min_mass_index]
    max_mass = mass[max_mass_index]
    den = (star.mass - min_mass) / (max_mass - min_mass)

    if case_1 == 0 and case_2 == 0:
        ym1 = y1 + (y3 - y1) * den
        ym2 = y2 + (y4 - y2) * den
        xm1 = x1 + (x3 - x1) * den
        xm2 = x2 + (x4 - x2) * den
        return xm1 + (star.cooling_time - ym1) / (ym2 - ym1) * (xm2 - xm1)

    if case_1 == 0 and case_2 == 1:
        xm1 = x1 + (x2 - x1) * (star.cooling_time - y1) / (y2 - y1)
        return xm1 + (x3 - xm1) * den

    if case_1 == 1 and case_2 == 0:
        xm2 = x3 + (x4 - x3) * (star.cooling_time - y3) / (y4 - y3)
        return x1 + (xm2 - x1) * den

    return x1 + (x3 - x1) * den


def get_xm(star: Star,
           cooling_time: np.ndarray,
           pre_wd_lifetime: np.ndarray,
           interest_sequence: np.ndarray,
           rows_counts: np.ndarray,
           mass_index: int,
           by_logarithm: bool,
           one_model: bool = False) -> float:
    if star.cooling_time < cooling_time[mass_index, 0]:
        return get_extrapolated_xm(
            star=star,
            cooling_time=cooling_time,
            pre_wd_lifetime=pre_wd_lifetime,
            interest_sequence=interest_sequence,
            mass_index=mass_index,
            min_row_index=0,
            by_logarithm=by_logarithm,
            one_model=one_model)

    rows_count = rows_counts[mass_index]

    if star.cooling_time > cooling_time[mass_index, rows_count]:
        return get_extrapolated_xm(
            star=star,
            cooling_time=cooling_time,
            pre_wd_lifetime=pre_wd_lifetime,
            interest_sequence=interest_sequence,
            mass_index=mass_index,
            min_row_index=rows_count - 1,
            by_logarithm=by_logarithm,
            one_model=one_model)

    for row_index in range(rows_count - 1):
        if (cooling_time[mass_index, row_index] <= star.cooling_time
                <= cooling_time[mass_index, row_index + 1]):
            y1 = cooling_time[mass_index, row_index]
            y2 = cooling_time[mass_index, row_index + 1]
            x1 = interest_sequence[mass_index, row_index]
            x2 = interest_sequence[mass_index, row_index + 1]
            deltf = (star.cooling_time - y1) / (y2 - y1)

            return x1 + deltf * (x2 - x1)


def get_extrapolated_xm(star: Star,
                        cooling_time: np.ndarray,
                        pre_wd_lifetime: np.ndarray,
                        interest_sequence: np.ndarray,
                        mass_index: int,
                        min_row_index: int,
                        by_logarithm: bool,
                        one_model: bool = False) -> float:
    if one_model:
        deltf = (cooling_time[mass_index, min_row_index + 1]
                 - cooling_time[mass_index, min_row_index])
        s = ((interest_sequence[mass_index, min_row_index + 1]
              - interest_sequence[mass_index, min_row_index]) / deltf)
        b = (interest_sequence[mass_index, min_row_index + 1]
             - s * cooling_time[mass_index, min_row_index + 1])
        return s * star.cooling_time + b

    deltf = log10((cooling_time[mass_index, min_row_index + 1]
                   + pre_wd_lifetime[mass_index])
                  / (cooling_time[mass_index, min_row_index]
                     + pre_wd_lifetime[mass_index]))

    if by_logarithm:
        s = log10(interest_sequence[mass_index, min_row_index + 1]
                  / interest_sequence[mass_index, min_row_index]) / deltf
        b = (log10(interest_sequence[mass_index, min_row_index + 1])
             - s * log10(cooling_time[mass_index, min_row_index + 1]
                         + pre_wd_lifetime[mass_index]))
        return 10.0 ** (s * log10(star.cooling_time
                                  + pre_wd_lifetime[mass_index]) + b)

    s = ((interest_sequence[mass_index, min_row_index + 1]
          - interest_sequence[mass_index, min_row_index]) / deltf)
    b = (interest_sequence[mass_index, min_row_index + 1]
         - s * log10(cooling_time[mass_index, min_row_index + 1]
                     + pre_wd_lifetime[mass_index]))
    return s * log10(star.cooling_time + pre_wd_lifetime[mass_index]) + b


def interpolate_magnitudes(star: Star,
                           color_table: Dict[str, np.ndarray],
                           luminosity: float) -> Tuple[float, ...]:
    if star.mass <= color_table['mass'][0]:
        rows_count_1 = color_table['rows_counts'][0]
        rows_count_2 = color_table['rows_counts'][1]
        min_mass_index = 0
    elif star.mass > color_table['mass'][-1]:
        rows_count_1 = color_table['rows_counts'][-2]
        rows_count_2 = color_table['rows_counts'][-1]
        min_mass_index = color_table['mass'].size() - 1
    else:
        for mass_index in range(color_table['mass'].size() - 1):
            if (color_table['mass'][mass_index] < star.mass
                    <= color_table['mass'][mass_index + 1]):
                rows_count_1 = color_table['rows_counts'][mass_index]
                rows_count_2 = color_table['rows_counts'][mass_index + 1]
                min_mass_index = mass_index
                break

    return get_magnitudes(star=star,
                          luminosity=luminosity,
                          color_table=color_table,
                          rows_count_1=rows_count_1,
                          rows_count_2=rows_count_2,
                          min_mass_index=min_mass_index)


def get_magnitudes(star: Star,
                   luminosity: float,
                   color_table: Dict[str, np.ndarray],
                   rows_count_1: int,
                   rows_count_2: int,
                   min_mass_index: int
                   ) -> Tuple[float, ...]:
    if (luminosity > color_table['luminosity'][min_mass_index, 0]
            or luminosity > color_table['luminosity'][min_mass_index + 1, 0]):
        return get_extrapolated_magnitudes_by_luminosity(
            star=star,
            luminosity=luminosity,
            color_table=color_table,
            row_index_1=0,
            row_index_2=0,
            mass_index=min_mass_index)

    if (luminosity < color_table['luminosity'][min_mass_index, rows_count_1]
        or luminosity < color_table['luminosity'][min_mass_index + 1,
                                                  rows_count_2]):
        return get_extrapolated_magnitudes_by_luminosity(
            star=star,
            luminosity=luminosity,
            color_table=color_table,
            row_index_1=rows_count_1,
            row_index_2=rows_count_2,
            mass_index=min_mass_index)

    return get_interpolated_magnitudes_by_luminosity(
        star=star,
        rows_count_1=rows_count_1,
        rows_count_2=rows_count_2,
        color_table=color_table,
        luminosity=luminosity,
        mass_index=min_mass_index)


def get_interpolated_magnitudes_by_luminosity(
        star: Star,
        rows_count_1: int,
        rows_count_2: int,
        color_table: Dict[str, np.ndarray],
        luminosity: float,
        mass_index: int) -> Tuple[float, ...]:
    check_1 = 0
    check_2 = 0

    for row_index in range(rows_count_1 - 1):
        if (color_table['luminosity'][mass_index, row_index + 1] <= luminosity
                <= color_table['luminosity'][mass_index, row_index]):
            row_index_1 = row_index
            check_1 = 1
            break

    for row_index in range(rows_count_2 - 1):
        if (color_table['luminosity'][mass_index + 1, row_index + 1]
                <= luminosity
                <= color_table['luminosity'][mass_index + 1, row_index]):
            row_index_2 = row_index
            check_2 = 1
            break

    if check_1 == 1 and check_2 == 1:
        a1 = luminosity - color_table['luminosity'][mass_index, row_index_1]
        a2 = luminosity - color_table['luminosity'][mass_index + 1,
                                                    row_index_2]
        b1 = (color_table['luminosity'][mass_index, row_index_1 + 1]
              - color_table['luminosity'][mass_index, row_index_1])
        b2 = (color_table['luminosity'][mass_index + 1, row_index_2 + 1]
              - color_table['luminosity'][mass_index + 1, row_index_2])
        u_ubvri_absolute = get_interpolated_magnitude(
            star=star,
            table_magnitude=color_table['u_ubvri_absolute'],
            table_mass=color_table['mass'],
            row_index_1=row_index_1,
            row_index_2=row_index_2,
            mass_index=mass_index,
            a1=a1,
            a2=a2,
            b1=b1,
            b2=b2)
        b_ubvri_absolute = get_interpolated_magnitude(
            star=star,
            table_magnitude=color_table['b_ubvri_absolute'],
            table_mass=color_table['mass'],
            row_index_1=row_index_1,
            row_index_2=row_index_2,
            mass_index=mass_index,
            a1=a1,
            a2=a2,
            b1=b1,
            b2=b2)
        v_ubvri_absolute = get_interpolated_magnitude(
            star=star,
            table_magnitude=color_table['v_ubvri_absolute'],
            table_mass=color_table['mass'],
            row_index_1=row_index_1,
            row_index_2=row_index_2,
            mass_index=mass_index,
            a1=a1,
            a2=a2,
            b1=b1,
            b2=b2)
        r_ubvri_absolute = get_interpolated_magnitude(
            star=star,
            table_magnitude=color_table['r_ubvri_absolute'],
            table_mass=color_table['mass'],
            row_index_1=row_index_1,
            row_index_2=row_index_2,
            mass_index=mass_index,
            a1=a1,
            a2=a2,
            b1=b1,
            b2=b2)
        i_ubvri_absolute = get_interpolated_magnitude(
            star=star,
            table_magnitude=color_table['i_ubvri_absolute'],
            table_mass=color_table['mass'],
            row_index_1=row_index_1,
            row_index_2=row_index_2,
            mass_index=mass_index,
            a1=a1,
            a2=a2,
            b1=b1,
            b2=b2)

        return (u_ubvri_absolute,
                b_ubvri_absolute,
                v_ubvri_absolute,
                r_ubvri_absolute,
                i_ubvri_absolute)


def get_extrapolated_magnitudes_by_luminosity(
        star: Star,
        luminosity: float,
        color_table: Dict[str, np.ndarray],
        row_index_1: int,
        row_index_2: int,
        mass_index: int
        ) -> Tuple[float, ...]:
    u_ubvri_absolute = get_abs_magnitude(
        star=star,
        luminosity=luminosity,
        table_absolute_magnitude=color_table['u_ubvri_absolute'],
        table_luminosity=color_table['luminosity'],
        table_mass=color_table['mass'],
        row_index_1=row_index_1,
        row_index_2=row_index_2,
        mass_index=mass_index)
    b_ubvri_absolute = get_abs_magnitude(
        star=star,
        luminosity=luminosity,
        table_absolute_magnitude=color_table['b_ubvri_absolute'],
        table_luminosity=color_table['luminosity'],
        table_mass=color_table['mass'],
        row_index_1=row_index_1,
        row_index_2=row_index_2,
        mass_index=mass_index)
    v_ubvri_absolute = get_abs_magnitude(
        star=star,
        luminosity=luminosity,
        table_absolute_magnitude=color_table['v_ubvri_absolute'],
        table_luminosity=color_table['luminosity'],
        table_mass=color_table['mass'],
        row_index_1=row_index_1,
        row_index_2=row_index_2,
        mass_index=mass_index)
    r_ubvri_absolute = get_abs_magnitude(
        star=star,
        luminosity=luminosity,
        table_absolute_magnitude=color_table['r_ubvri_absolute'],
        table_luminosity=color_table['luminosity'],
        table_mass=color_table['mass'],
        row_index_1=row_index_1,
        row_index_2=row_index_2,
        mass_index=mass_index)
    i_ubvri_absolute = get_abs_magnitude(
        star=star,
        luminosity=luminosity,
        table_absolute_magnitude=color_table['i_ubvri_absolute'],
        table_luminosity=color_table['luminosity'],
        table_mass=color_table['mass'],
        row_index_1=row_index_1,
        row_index_2=row_index_2,
        mass_index=mass_index)

    return (u_ubvri_absolute,
            b_ubvri_absolute,
            v_ubvri_absolute,
            r_ubvri_absolute,
            i_ubvri_absolute)


def get_interpolated_magnitude(star: Star,
                               table_magnitude: np.ndarray,
                               table_mass: np.ndarray,
                               row_index_1: int,
                               row_index_2: int,
                               mass_index: int,
                               a1: float,
                               a2: float,
                               b1: float,
                               b2: float) -> float:
    c_1 = (table_magnitude[mass_index, row_index_1]
           + (table_magnitude[mass_index, row_index_1 + 1]
              - table_magnitude[mass_index, row_index_1]) * a1
           / b1)
    c_2 = (table_magnitude[mass_index + 1, row_index_2]
           + (table_magnitude[mass_index + 1, row_index_2 + 1]
              - table_magnitude[mass_index + 1, row_index_2]) * a2
           / b2)
    return extrapolate(luminosity_or_mass=star.mass,
                       interest_value_1=c_1,
                       interest_value_2=c_2,
                       table_luminosity_or_mass_1=table_mass[0],
                       table_luminosity_or_mass_2=table_mass[1])


def get_abs_magnitude(star: Star,
                      luminosity: float,
                      table_absolute_magnitude: np.ndarray,
                      table_luminosity: np.ndarray,
                      table_mass: np.ndarray,
                      row_index_1: int,
                      row_index_2: int,
                      mass_index: int) -> float:
    c_1 = extrapolate(
        luminosity_or_mass=luminosity,
        interest_value_1=table_absolute_magnitude[mass_index, row_index_1],
        interest_value_2=table_absolute_magnitude[mass_index, row_index_1 + 1],
        table_luminosity_or_mass_1=table_luminosity[mass_index, row_index_1],
        table_luminosity_or_mass_2=table_luminosity[mass_index,
                                                    row_index_1 + 1])
    c_2 = extrapolate(
        luminosity_or_mass=luminosity,
        interest_value_1=table_absolute_magnitude[mass_index + 1,
                                                  row_index_2],
        interest_value_2=table_absolute_magnitude[mass_index + 1,
                                                  row_index_2 + 1],
        table_luminosity_or_mass_1=table_luminosity[mass_index + 1,
                                                    row_index_2],
        table_luminosity_or_mass_2=table_luminosity[mass_index + 1,
                                                    row_index_2 + 1])
    abs_magnitude = extrapolate(
        luminosity_or_mass=star.mass,
        interest_value_1=c_1,
        interest_value_2=c_2,
        table_luminosity_or_mass_1=table_mass[mass_index],
        table_luminosity_or_mass_2=table_mass[mass_index + 1])
    if abs_magnitude < 0.:
        abs_magnitude = 0.

    return abs_magnitude


def extrapolate(luminosity_or_mass: float,
                interest_value_1: float,
                interest_value_2: float,
                table_luminosity_or_mass_1: float,
                table_luminosity_or_mass_2: float) -> float:
    s = (interest_value_2 - interest_value_1) / (table_luminosity_or_mass_2
                                                 - table_luminosity_or_mass_1)
    b = interest_value_1 - s * table_luminosity_or_mass_1
    return s * luminosity_or_mass + b

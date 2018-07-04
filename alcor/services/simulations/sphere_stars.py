import pandas as pd
import sympy as sym

from alcor.models.star import GalacticDiskType
from alcor.services import inverse_transform


# TODO: take these consts to settings file
def generate_stars(*,
                   thin_disk_age: float,
                   thick_disk_age: float,
                   halo_age: float,
                   halo_stars_formation_time: float,
                   burst_age: float,
                   thick_disk_formation_rate_exponent: float,
                   thick_disk_stars_fraction: float,
                   halo_stars_fraction: float,
                   initial_mass_function_exponent: float,
                   burst_formation_factor: float,
                   thin_disk_stars_count: int,
                   min_mass: float,
                   max_mass: float) -> pd.DataFrame:
    max_age = max(thin_disk_age, thick_disk_age, halo_age)

    thin_disk_birth_times = thin_disk_birth_time(
            initial_time=max_age - thin_disk_age,
            burst_initial_time=max_age - burst_age,
            max_age=max_age,
            burst_factor=burst_formation_factor)
    birth_times = inverse_transform.sample(thin_disk_birth_times,
                                           size=thin_disk_stars_count)
    thin_disk_stars = pd.DataFrame(dict(
            birth_time=birth_times,
            galactic_disk_type=GalacticDiskType.thin))

    shift = max_age - thick_disk_age
    thick_disk_birth_times = thick_disk_birth_time(
            shift=shift,
            max_age=max_age,
            exponent=thick_disk_formation_rate_exponent)
    thin_disk_stars_fraction = (1. - thick_disk_stars_fraction
                                - halo_stars_fraction)
    thick_disk_stars_count = int(thin_disk_stars.shape[0]
                                 * thick_disk_stars_fraction
                                 / thin_disk_stars_fraction)
    birth_times = inverse_transform.sample(thick_disk_birth_times,
                                           size=thick_disk_stars_count)
    thick_disk_stars = pd.DataFrame(dict(
            birth_time=birth_times,
            galactic_disk_type=GalacticDiskType.thick))

    birth_initial_time = max_age - halo_age
    birth_final_time = birth_initial_time + halo_stars_formation_time
    halo_birth_times = halo_birth_time(initial_time=birth_initial_time,
                                       final_time=birth_final_time)
    halo_stars_count = int(thin_disk_stars.shape[0] * halo_stars_fraction
                           / thin_disk_stars_fraction)
    birth_times = inverse_transform.sample(halo_birth_times,
                                           size=halo_stars_count)
    halo_stars = pd.DataFrame(dict(birth_time=birth_times,
                                   galactic_disk_type=GalacticDiskType.halo))

    stars = pd.concat([thin_disk_stars, thick_disk_stars, halo_stars])

    progenitor_masses = progenitor_mass(
            min_mass=min_mass,
            max_mass=max_mass,
            exponent=initial_mass_function_exponent)
    progenitors_masses = inverse_transform.sample(progenitor_masses,
                                                  size=stars.shape[0])
    stars['progenitor_mass'] = progenitors_masses

    return stars


def generate_halo_stars(*,
                        stars_count: int,
                        formation_time: float,
                        max_age: float,
                        halo_age: float) -> pd.DataFrame:
    x = sym.Symbol('x')
    birth_initial_time = max_age - halo_age
    birth_final_time = birth_initial_time + formation_time
    f = sym.Piecewise((0, x < birth_initial_time),
                      (1, x <= birth_final_time),
                      (0, True))
    birth_times = inverse_transform.sample(f, size=stars_count)

    return pd.DataFrame(dict(birth_time=birth_times,
                             galactic_disk_type=GalacticDiskType.halo))


def thin_disk_birth_time(*,
                         initial_time: float,
                         burst_initial_time: float,
                         max_age: float,
                         burst_factor: float) -> sym.Piecewise:
    x = sym.Symbol('x')
    return sym.Piecewise((0, x < initial_time),
                         (1, x < burst_initial_time),
                         (burst_factor, x <= max_age),
                         (0, True))


def thick_disk_birth_time(*,
                          shift: float,
                          max_age: float,
                          exponent: float) -> sym.Piecewise:
    x = sym.Symbol('x')
    return sym.Piecewise(
            (0., x <= shift),
            ((x - shift) * sym.exp(-(x - shift) / exponent), x <= max_age),
            (0., True))


def halo_birth_time(*,
                    initial_time: float,
                    final_time: float) -> sym.Piecewise:
    x = sym.Symbol('x')
    return sym.Piecewise((0, x < initial_time),
                         (1, x <= final_time),
                         (0, True))


def progenitor_mass(*,
                    min_mass: float,
                    max_mass: float,
                    exponent: float) -> sym.Piecewise:
    x = sym.Symbol('x')
    return sym.Piecewise((0, x < min_mass),
                         (x ** exponent, x < max_mass),
                         (0, True))

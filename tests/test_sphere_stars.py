import math
import random

import numpy as np
import pandas as pd

from alcor.services.simulations.sphere_stars import (
    halo_star_birth_time,
    thin_disk_star_birth_time,
    thick_disk_star_birth_time,
    initial_star_mass_by_salpeter,
    normalization_const,
    get_birth_rates,
    generate_thin_disk_stars,
    generate_thick_disk_stars,
    generate_halo_stars)

UNIVERSE_AGE = 14.


# TODO: find out how to test RNGs
def test_halo_star_birth_time(halo_birth_initial_time: float,
                              halo_stars_formation_time: float) -> None:
    birth_time = halo_star_birth_time(
            birth_initial_time=halo_birth_initial_time,
            formation_time=halo_stars_formation_time,
            generator=random.uniform)

    assert isinstance(birth_time, float)
    assert birth_time >= halo_birth_initial_time
    assert 0. <= birth_time <= birth_time + halo_stars_formation_time


def test_thin_disk_star_birth_time(bin_initial_time: float,
                                   time_increment: float) -> None:
    birth_time = thin_disk_star_birth_time(
            bin_initial_time=bin_initial_time,
            time_increment=time_increment,
            generator=random.uniform)

    assert isinstance(birth_time, float)
    assert UNIVERSE_AGE >= birth_time >= bin_initial_time


def test_thick_disk_star_birth_time(age: float,
                                    formation_rate_parameter: float,
                                    max_formation_rate: float,
                                    thick_disk_birth_initial_time: float
                                    ) -> None:
    birth_time = thick_disk_star_birth_time(
            age=age,
            formation_rate_parameter=formation_rate_parameter,
            max_formation_rate=max_formation_rate,
            birth_initial_time=thick_disk_birth_initial_time,
            generator=random.uniform    )

    assert isinstance(birth_time, float)
    assert UNIVERSE_AGE >= birth_time >= thick_disk_birth_initial_time


def test_initial_star_mass_by_salpeter(exponent: float,
                                       min_mass: float,
                                       max_mass: float) -> None:
    mass = initial_star_mass_by_salpeter(exponent=exponent,
                                         min_mass=min_mass,
                                         max_mass=max_mass,
                                         generator=random.uniform)

    assert isinstance(mass, float)
    assert min_mass <= mass <= max_mass


def test_normalization_const(star_formation_rate_param: float,
                             thin_disk_age_gyr: float,
                             sigma: float) -> None:
    normalization_constant = normalization_const(
            star_formation_rate_param=star_formation_rate_param,
            thin_disk_age_gyr=thin_disk_age_gyr,
            sigma=sigma)

    assert isinstance(normalization_constant, float)
    assert math.isfinite(normalization_constant)


def test_get_birth_rates(times: np.ndarray,
                         burst_init_time: float,
                         birth_rate: float,
                         burst_birth_rate: float) -> None:
    birth_rates = get_birth_rates(times=times,
                                  burst_init_time=burst_init_time,
                                  birth_rate=birth_rate,
                                  burst_birth_rate=burst_birth_rate)

    assert isinstance(birth_rates, np.ndarray)
    assert birth_rates.size == times.size


def test_generate_thin_disk_stars(max_age: float,
                                  time_bins_count: int,
                                  burst_age: float,
                                  initial_mass_function_parameter: float,
                                  thin_disk_age_gyr: float,
                                  max_stars_count: int,
                                  sector_radius_kpc: float,
                                  burst_formation_factor: float,
                                  star_formation_rate_param: float,
                                  mass_reduction_factor: float) -> None:
    thin_disk_stars = generate_thin_disk_stars(
            max_age=max_age,
            time_bins_count=time_bins_count,
            burst_age=burst_age,
            initial_mass_function_parameter=initial_mass_function_parameter,
            age=thin_disk_age_gyr,
            max_stars_count=max_stars_count,
            sector_radius_kpc=sector_radius_kpc,
            burst_formation_factor=burst_formation_factor,
            star_formation_rate_param=star_formation_rate_param,
            mass_reduction_factor=mass_reduction_factor,
            generator=random.uniform)

    assert isinstance(thin_disk_stars, pd.DataFrame)
    assert thin_disk_stars.columns.size > 0


def test_generate_thick_disk_stars(thin_disk_stars_fraction: float,
                                   thick_disk_stars_fraction: float,
                                   thin_disk_stars_count: int,
                                   initial_mass_function_parameter: float,
                                   thick_disk_age: float,
                                   max_age: float,
                                   thick_disk_sfr_param) -> None:
    thick_disk_stars = generate_thick_disk_stars(
            thin_disk_stars_fraction=thin_disk_stars_fraction,
            thick_disk_stars_fraction=thick_disk_stars_fraction,
            thin_disk_stars_count=thin_disk_stars_count,
            initial_mass_function_parameter=initial_mass_function_parameter,
            age=thick_disk_age,
            max_age=max_age,
            formation_rate_parameter=thick_disk_sfr_param,
            generator=random.uniform)

    assert isinstance(thick_disk_stars, pd.DataFrame)
    assert thick_disk_stars.columns.size > 0


def test_generate_halo_stars(thin_disk_stars_count: int,
                             halo_stars_fraction: float,
                             thin_disk_stars_fraction: float,
                             initial_mass_function_parameter: float,
                             max_age: float,
                             halo_age: float,
                             halo_stars_formation_time: float) -> None:
    halo_stars = generate_halo_stars(
            thin_disk_stars_count=thin_disk_stars_count,
            halo_stars_fraction=halo_stars_fraction,
            thin_disk_stars_fraction=thin_disk_stars_fraction,
            initial_mass_function_parameter=initial_mass_function_parameter,
            max_age=max_age,
            halo_age=halo_age,
            formation_time=halo_stars_formation_time,
            generator=random.uniform)

    assert isinstance(halo_stars, pd.DataFrame)
    assert halo_stars.columns.size > 0

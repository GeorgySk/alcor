import math
import random
from typing import Callable

import numpy as np
import pandas as pd

from alcor.services.simulations.sphere_stars import (
    halo_star_birth_time,
    thin_disk_star_birth_time,
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


def test_normalization_const(formation_rate_parameter: float,
                             disk_age: float,
                             sigma: float) -> None:
    normalization_constant = normalization_const(
            formation_rate_parameter=formation_rate_parameter,
            thin_disk_age=disk_age,
            sigma=sigma)

    assert isinstance(normalization_constant, float)
    assert math.isfinite(normalization_constant)


def test_get_birth_rates(times: np.ndarray,
                         burst_initial_time: float,
                         birth_rate: float,
                         burst_birth_rate: float) -> None:
    birth_rates = get_birth_rates(times=times,
                                  burst_initial_time=burst_initial_time,
                                  birth_rate=birth_rate,
                                  burst_birth_rate=burst_birth_rate)

    assert isinstance(birth_rates, np.ndarray)
    assert birth_rates.size == times.size


def test_generate_thin_disk_stars(
        max_age: float,
        time_bins_count: int,
        burst_age: float,
        disk_age: float,
        max_stars_count: int,
        sector_radius_kpc: float,
        burst_formation_factor: float,
        formation_rate_parameter: float,
        mass_reduction_factor: float,
        initial_mass_generator: Callable[[float], float]) -> None:
    thin_disk_stars = generate_thin_disk_stars(
            max_age=max_age,
            time_bins_count=time_bins_count,
            burst_age=burst_age,
            disk_age=disk_age,
            max_stars_count=max_stars_count,
            sector_radius_kpc=sector_radius_kpc,
            burst_formation_factor=burst_formation_factor,
            formation_rate_parameter=formation_rate_parameter,
            mass_reduction_factor=mass_reduction_factor,
            generator=random.uniform,
            initial_mass_generator=initial_mass_generator)

    assert isinstance(thin_disk_stars, pd.DataFrame)
    assert thin_disk_stars.columns.size > 0


def test_generate_thick_disk_stars(
        stars_count: int,
        thick_disk_age: float,
        formation_rate_exponent: float,
        birth_initial_time: float,
        initial_mass_generator: Callable[[float], float]) -> None:
    thick_disk_stars = generate_thick_disk_stars(
            stars_count=stars_count,
            disk_age=thick_disk_age,
            birth_initial_time=birth_initial_time,
            formation_rate_exponent=formation_rate_exponent,
            generator=random.uniform,
            initial_mass_generator=initial_mass_generator)

    assert isinstance(thick_disk_stars, pd.DataFrame)
    assert thick_disk_stars.columns.size > 0


def test_generate_halo_stars(
        stars_count: int,
        birth_initial_time: float,
        halo_stars_formation_time: float,
        initial_mass_generator: Callable[[float], float]) -> None:
    halo_stars = generate_halo_stars(
            stars_count=stars_count,
            birth_initial_time=birth_initial_time,
            formation_time=halo_stars_formation_time,
            generator=random.uniform,
            initial_mass_generator=initial_mass_generator)

    assert isinstance(halo_stars, pd.DataFrame)
    assert halo_stars.columns.size > 0

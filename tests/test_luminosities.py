import math
import numpy as np
import pandas as pd

from alcor.services.simulations.luminosities import (
    set_metallicities,
    extrapolated_times,
    estimated_times,
    estimate_lifetime,
    get_main_sequence_lifetimes,
    set_cooling_times,
    set_masses,
    white_dwarf_masses,
    get_white_dwarfs)


def test_set_metallicities(stars_w_galactic_disk_types: pd.DataFrame,
                           solar_metallicity: float,
                           subsolar_metallicity: float) -> None:
    set_metallicities(stars_w_galactic_disk_types,
                      subsolar_metallicity=subsolar_metallicity,
                      solar_metallicity=solar_metallicity)

    thin_disk_stars_mask = (stars_w_galactic_disk_types['galactic_disk_type']
                            == 'thin')
    thin_disk_stars = stars_w_galactic_disk_types[thin_disk_stars_mask]
    thin_disk_stars_metallicities = thin_disk_stars['metallicity']
    thick_disk_stars_mask = (stars_w_galactic_disk_types['galactic_disk_type']
                             == 'thick')
    thick_disk_stars = stars_w_galactic_disk_types[thick_disk_stars_mask]
    thick_disk_stars_metallicities = thick_disk_stars['metallicity']
    halo_stars_mask = (stars_w_galactic_disk_types['galactic_disk_type']
                       == 'halo')
    halo_stars = stars_w_galactic_disk_types[halo_stars_mask]
    halo_stars_metallicities = halo_stars['metallicity']

    assert isinstance(stars_w_galactic_disk_types, pd.DataFrame)
    assert stars_w_galactic_disk_types.shape[0] > 0
    assert 'galactic_disk_type' in stars_w_galactic_disk_types.columns
    assert 'metallicity' in stars_w_galactic_disk_types.columns
    assert (thin_disk_stars_metallicities == solar_metallicity).all()
    assert (thick_disk_stars_metallicities == solar_metallicity).all()
    assert (halo_stars_metallicities == subsolar_metallicity).all()


def test_extrapolated_times(masses: np.ndarray,
                            rightmost_mass: float,
                            rightmost_time: float) -> None:
    times = extrapolated_times(masses=masses,
                               rightmost_mass=rightmost_mass,
                               rightmost_time=rightmost_time)

    assert isinstance(times, np.ndarray)
    assert masses.shape == times.shape
    assert (times > 0.).all()


def test_estimated_times(masses: np.ndarray) -> None:
    model_solar_masses = np.array([1.00, 1.50, 1.75, 2.00, 2.25,
                                   2.50, 3.00, 3.50, 4.00, 5.00])
    model_solar_times = np.array([8.614, 1.968, 1.249, 0.865, 0.632,
                                  0.480, 0.302, 0.226, 0.149, 0.088])
    model_subsolar_masses = np.array([0.85, 1.00, 1.25, 1.50,
                                      1.75, 2.00, 3.00])
    model_subsolar_times = np.array([10.34, 5.756, 2.623, 1.412,
                                     0.905, 0.639, 0.245])

    solar_times = estimated_times(masses=masses,
                                  model_masses=model_solar_masses,
                                  model_times=model_solar_times)
    subsolar_times = estimated_times(masses=masses,
                                     model_masses=model_subsolar_masses,
                                     model_times=model_subsolar_times)

    assert isinstance(solar_times, np.ndarray)
    assert isinstance(subsolar_times, np.ndarray)
    assert masses.shape == solar_times.shape
    assert masses.shape == subsolar_times.shape


def test_estimate_lifetime() -> None:
    estimated_lifetime = estimate_lifetime(metallicity=3.,
                                           subsolar_main_sequence_lifetime=1.,
                                           solar_main_sequence_lifetime=2.,
                                           subsolar_metallicity=1.,
                                           solar_metallicity=2.)

    assert math.isclose(estimated_lifetime, 3.)


def test_get_main_sequence_lifetimes(masses: np.ndarray,
                                     metallicities: np.ndarray) -> None:
    main_sequence_lifetimes = get_main_sequence_lifetimes(
            masses=masses,
            metallicities=metallicities,
            solar_metallicity=0.01,
            subsolar_metallicity=0.001)

    assert isinstance(main_sequence_lifetimes, np.ndarray)


def test_set_cooling_times(stars_without_cooling_times: pd.DataFrame) -> None:
    set_cooling_times(stars_without_cooling_times,
                      max_galactic_structure_age=12.,
                      subsolar_metallicity=0.001,
                      solar_metallicity=0.01)

    assert isinstance(stars_without_cooling_times, pd.DataFrame)
    assert 'cooling_time' in stars_without_cooling_times.columns


def test_get_white_dwarf_masses(progenitor_masses: np.ndarray) -> None:
    masses = white_dwarf_masses(progenitor_masses)

    assert isinstance(masses, np.ndarray)
    assert progenitor_masses.shape[0] == masses.shape[0]


def test_set_masses(stars_without_masses: pd.DataFrame) -> None:
    columns_before = stars_without_masses.columns
    set_masses(stars_without_masses,
               ifmr_parameter=1.)
    columns_after = stars_without_masses.columns

    assert isinstance(stars_without_masses, pd.DataFrame)
    assert len(columns_after) == len(columns_before) + 1
    assert 'mass' not in columns_before
    assert 'mass' in columns_after


def test_get_white_dwarfs(main_sequence_stars: pd.DataFrame) -> None:
    white_dwarfs = get_white_dwarfs(main_sequence_stars,
                                    max_galactic_structure_age=12.)

    assert isinstance(white_dwarfs, pd.DataFrame)
    assert 'metallicity' in white_dwarfs.columns
    assert 'cooling_time' in white_dwarfs.columns
    assert 'mass' in white_dwarfs.columns

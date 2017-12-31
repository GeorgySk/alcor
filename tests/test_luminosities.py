import math
import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

from alcor.services.simulations import luminosities


# TODO: add tests checking that we didn't mutate anything
from alcor.services.simulations.luminosities import get_metallicities


def test_extrapolated_times(masses: np.ndarray,
                            rightmost_mass: float,
                            rightmost_time: float) -> None:
    times = luminosities.extrapolate_main_sequence_lifetimes(
            masses=masses,
            rightmost_mass=rightmost_mass,
            rightmost_time=rightmost_time)

    assert isinstance(times, np.ndarray)
    assert masses.shape == times.shape
    assert (times > 0.).all()


def test_estimated_times(masses: np.ndarray,
                         rightmost_mass: float,
                         rightmost_time: float,
                         spline: InterpolatedUnivariateSpline) -> None:
    masses_before = np.copy(masses)
    solar_times = luminosities.estimated_times(
            masses=masses,
            rightmost_model_mass=rightmost_mass,
            rightmost_model_time=rightmost_time,
            spline=spline)

    assert isinstance(solar_times, np.ndarray)
    assert masses.shape == solar_times.shape
    assert (masses_before == masses).all()


def test_estimate_lifetime() -> None:
    estimated_lifetime = luminosities.estimate_lifetime(
            metallicity=3.,
            subsolar_main_sequence_lifetime=1.,
            solar_main_sequence_lifetime=2.,
            subsolar_metallicity=1.,
            solar_metallicity=2.)

    assert math.isclose(estimated_lifetime, 3.)


def test_get_main_sequence_lifetimes(masses: np.ndarray,
                                     metallicities: np.ndarray,
                                     solar_metallicity: float,
                                     subsolar_metallicity: float) -> None:
    main_sequence_lifetimes = luminosities.main_sequence_stars_lifetimes(
            masses=masses,
            metallicities=metallicities,
            solar_metallicity=solar_metallicity,
            subsolar_metallicity=subsolar_metallicity)

    assert isinstance(main_sequence_lifetimes, np.ndarray)


def test_get_white_dwarf_masses(progenitor_masses: np.ndarray) -> None:
    progenitor_masses_before = np.copy(progenitor_masses)
    masses = luminosities.white_dwarf_masses(progenitor_masses)

    assert isinstance(masses, np.ndarray)
    assert progenitor_masses.shape[0] == masses.shape[0]
    assert (progenitor_masses_before == progenitor_masses).all()


def test_get_white_dwarfs(main_sequence_stars: pd.DataFrame) -> None:
    # TODO: after merge there will be a fixture for ages
    white_dwarfs = luminosities.white_dwarfs(main_sequence_stars,
                                             max_galactic_structure_age=12.)

    assert isinstance(white_dwarfs, pd.DataFrame)
    assert 'metallicity' in white_dwarfs.columns
    assert 'cooling_time' in white_dwarfs.columns
    assert 'mass' in white_dwarfs.columns


def test_get_metallicities(galactic_disks_types: np.ndarray,
                           subsolar_metallicity: float,
                           solar_metallicity: float) -> None:
    galactic_disks_types_before = np.copy(galactic_disks_types)
    metallicities = get_metallicities(
            galactic_disks_types=galactic_disks_types,
            subsolar_metallicity=subsolar_metallicity,
            solar_metallicity=solar_metallicity)

    assert isinstance(metallicities, np.ndarray)
    assert (galactic_disks_types_before == galactic_disks_types).all()
    assert metallicities.size == galactic_disks_types.size

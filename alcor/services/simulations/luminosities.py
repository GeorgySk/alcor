import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline


def get_white_dwarfs(stars: pd.DataFrame,
                     *,
                     max_galactic_structure_age: float,
                     ifmr_parameter: float = 1.,
                     chandrasekhar_limit: float = 1.4,
                     max_mass: float = 10.5,
                     solar_metallicity: float = 0.01,
                     subsolar_metallicity: float = 0.001) -> pd.DataFrame:
    """Filter white dwarfs stars(WD) from initial sample of main sequence stars
    and assign metallicities, cooling times and masses.

    Keyword arguments:
    stars -- data frame with main sequence stars
    max_galactic_structure_age -- the highest age of thin disk, thick disk
                                  and halo
    ifmr_parameter -- factor by which WD mass (calculated from Initial-to-Final
                      Mass Relation (IFMR)) is multiplied (default 1.)
    chandrasekhar_limit -- maximum mass of a stable WD (default 1.4)
    max_mass -- maximum mass of a main sequence star that can generate a WD
                (default 10.5)
    solar_metallicity -- metallicity assigned to all thin and thick disks WDs
                         due to relatively young ages (default 0.01)
    subsolar_metallicity -- metallicity assigned to all halo WDs
                            (default 0.001)
    """
    stars = filter_by_max_mass(stars,
                               max_mass=max_mass)
    set_metallicities(stars,
                      subsolar_metallicity=subsolar_metallicity,
                      solar_metallicity=solar_metallicity)
    set_cooling_times(stars,
                      max_galactic_structure_age=max_galactic_structure_age,
                      subsolar_metallicity=subsolar_metallicity,
                      solar_metallicity=solar_metallicity)
    stars = filter_by_cooling_time(stars)
    set_masses(stars,
               ifmr_parameter=ifmr_parameter)

    return filter_by_chandrasekhar_limit(stars,
                                         limit=chandrasekhar_limit)


def filter_by_max_mass(stars: pd.DataFrame,
                       *,
                       max_mass: float) -> pd.DataFrame:
    return stars[stars['progenitor_mass'] < max_mass]


def set_metallicities(stars: pd.DataFrame,
                      *,
                      subsolar_metallicity: float,
                      solar_metallicity: float) -> None:
    stars['metallicity'] = np.empty(stars.shape[0])

    halo_stars_mask = (stars['galactic_disk_type'] == 'halo')
    non_halo_stars_mask = (stars['galactic_disk_type'] != 'halo')

    stars.loc[halo_stars_mask, 'metallicity'] = subsolar_metallicity
    stars.loc[non_halo_stars_mask, 'metallicity'] = solar_metallicity


def set_cooling_times(stars: pd.DataFrame,
                      *,
                      max_galactic_structure_age: float,
                      solar_metallicity: float,
                      subsolar_metallicity: float) -> None:
    main_sequence_lifetimes = get_main_sequence_lifetimes(
            masses=stars['progenitor_mass'],
            metallicities=stars['metallicity'],
            solar_metallicity=solar_metallicity,
            subsolar_metallicity=subsolar_metallicity)

    stars['cooling_time'] = (max_galactic_structure_age - stars['birth_time']
                             - main_sequence_lifetimes)


def filter_by_cooling_time(stars: pd.DataFrame,
                           *,
                           min_cooling_time: float = 0.) -> pd.DataFrame:
    return stars[stars['cooling_time'] > min_cooling_time]


def set_masses(stars: pd.DataFrame,
               *,
               ifmr_parameter: float) -> None:
    stars['mass'] = ifmr_parameter * get_white_dwarf_masses(
            progenitor_masses=stars['progenitor_mass'])


def filter_by_chandrasekhar_limit(stars: pd.DataFrame,
                                  *,
                                  limit: float) -> pd.DataFrame:
    return stars[stars['mass'] <= limit]


def get_main_sequence_lifetimes(*,
                                masses: np.ndarray,
                                metallicities: np.ndarray,
                                solar_metallicity: float,
                                subsolar_metallicity: float) -> np.ndarray:
    """Calculate lifetime of a main sequence star
    according to model by Leandro & Renedo et al.(2010).
    Solar metallicity values from Althaus priv. comm (X = 0.725, Y = 0.265)
    Sub-solar metallicity values from Althaus priv. comm (X = 0.752, Y = 0.247)
    """
    model_solar_masses = np.array([1.00, 1.50, 1.75, 2.00, 2.25,
                                   2.50, 3.00, 3.50, 4.00, 5.00])
    model_solar_times = np.array([8.614, 1.968, 1.249, 0.865, 0.632,
                                  0.480, 0.302, 0.226, 0.149, 0.088])
    model_subsolar_masses = np.array([0.85, 1.00, 1.25, 1.50,
                                      1.75, 2.00, 3.00])
    model_subsolar_times = np.array([10.34, 5.756, 2.623, 1.412,
                                     0.905, 0.639, 0.245])

    solar_main_sequence_lifetimes = estimated_times(
            masses=masses,
            model_masses=model_solar_masses,
            model_times=model_solar_times)
    subsolar_main_sequence_lifetimes = estimated_times(
            masses=masses,
            model_masses=model_subsolar_masses,
            model_times=model_subsolar_times)

    estimate_lifetimes = np.vectorize(estimate_lifetime)

    return estimate_lifetimes(
            metallicity=metallicities,
            subsolar_main_sequence_lifetime=subsolar_main_sequence_lifetimes,
            solar_main_sequence_lifetime=solar_main_sequence_lifetimes,
            subsolar_metallicity=subsolar_metallicity,
            solar_metallicity=solar_metallicity)


def estimated_times(*,
                    masses: np.ndarray,
                    model_masses: np.ndarray,
                    model_times: np.ndarray) -> np.ndarray:
    times = np.empty(masses.shape[0])

    spline = InterpolatedUnivariateSpline(x=model_masses,
                                          y=model_times,
                                          k=1)

    masses_lt_max_mask = masses < model_masses[-1]
    masses_ge_max_mask = ~masses_lt_max_mask
    times[masses_lt_max_mask] = spline(masses[masses_lt_max_mask])
    times[masses_ge_max_mask] = extrapolated_times(
            masses=masses[masses_ge_max_mask],
            rightmost_mass=model_masses[-1],
            rightmost_time=model_times[-1])

    return times


def extrapolated_times(*,
                       masses: np.ndarray,
                       rightmost_mass: float,
                       rightmost_time: float) -> np.ndarray:
    """Extrapolate main sequence stars (progenitors) lifetime vs mass to the
    right. Make sure that no negative values will be produced."""
    return rightmost_time * rightmost_mass / masses


def estimate_lifetime(metallicity: float,
                      subsolar_main_sequence_lifetime: float,
                      solar_main_sequence_lifetime: float,
                      subsolar_metallicity: float,
                      solar_metallicity: float) -> float:
    spline = InterpolatedUnivariateSpline(
            x=(subsolar_metallicity, solar_metallicity),
            y=(subsolar_main_sequence_lifetime, solar_main_sequence_lifetime),
            k=1)
    return spline(metallicity)


def get_white_dwarf_masses(progenitor_masses: np.ndarray) -> np.ndarray:
    masses = np.empty(progenitor_masses.shape[0])

    low_progenitor_masses_mask = progenitor_masses < 2.7
    low_progenitor_masses = progenitor_masses[low_progenitor_masses_mask]

    medium_progenitor_masses_mask = ((progenitor_masses >= 2.7)
                                     & (progenitor_masses <= 6.))
    medium_progenitor_masses = progenitor_masses[medium_progenitor_masses_mask]

    high_progenitor_masses_mask = (progenitor_masses > 6.)
    high_progenitor_masses = progenitor_masses[high_progenitor_masses_mask]

    masses[low_progenitor_masses_mask] = 0.096 * low_progenitor_masses + 0.429
    masses[medium_progenitor_masses_mask] = (0.137 * medium_progenitor_masses
                                             + 0.3183)
    masses[high_progenitor_masses_mask] = (0.1057 * high_progenitor_masses
                                           + 0.5061)

    return masses

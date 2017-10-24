import numpy as np
import pandas as pd

from alcor.models.star import GalacticDiskType


def get_white_dwarfs(stars: pd.DataFrame,
                     max_galactic_structure_age: float,
                     ifmr_parameter: float,
                     chandrasekhar_limit: float = 1.4,
                     max_mass: float = 10.5,
                     solar_metallicity: float = 0.01,
                     subsolar_metallicity: float = 0.001) -> pd.DataFrame:
    valid_mass_stars = stars[stars['progenitor_mass'] < max_mass]

    # TODO: find out if this comparison will work
    halo_stars_mask = (valid_mass_stars['galactic_disk_type']
                       == GalacticDiskType.halo)
    halo_stars = valid_mass_stars[halo_stars_mask]
    non_halo_stars_mask = (valid_mass_stars['galactic_disk_type']
                           != GalacticDiskType.halo)
    non_halo_stars = valid_mass_stars[non_halo_stars_mask]

    halo_stars['metallicity'] = subsolar_metallicity
    non_halo_stars['metallicity'] = solar_metallicity

    stars = pd.concat(halo_stars, non_halo_stars)

    # TODO: fix this part
    main_sequence_lifetimes = []
    for progenitor_mass, metallicity in zip(stars['progenitor_mass'],
                                            stars['metallicity']):
        main_sequence_lifetimes.append(get_main_sequence_lifetime(
                mass=progenitor_mass,
                metallicity=metallicity))
    main_sequence_lifetimes = np.array(main_sequence_lifetimes)

    stars['cooling_times'] = (max_galactic_structure_age - stars['birth_time']
                              - main_sequence_lifetimes)

    valid_cooling_time_stars = stars[stars['cooling_times'] > 0.]

    valid_cooling_time_stars = ifmr_parameter * get_white_dwarf_masses(
            progenitor_masses=valid_cooling_time_stars['progenitor_mass'])

    white_dwarfs_mask = valid_cooling_time_stars['mass'] <= chandrasekhar_limit
    white_dwarfs = valid_cooling_time_stars[white_dwarfs_mask]

    return white_dwarfs


# TODO: use pandas
# According to model by Leandro & Renedo et al.(2010)
def get_main_sequence_lifetime(mass: float,
                               metallicity: float) -> float:
    main_sequence_masses = np.array([1.00, 1.50, 1.75, 2.00, 2.25,
                                     2.50, 3.00, 3.50, 4.00, 5.00])
    # Althaus priv. comm X = 0.725, Y = 0.265
    main_sequence_times = np.array([8.614, 1.968, 1.249, 0.865, 0.632,
                                    0.480, 0.302, 0.226, 0.149, 0.088])

    if mass < main_sequence_masses[0]:
        pen = ((main_sequence_times[1] - main_sequence_times[0])
               / (main_sequence_masses[1] - main_sequence_masses[0]))
        tsol = pen * mass + (main_sequence_times[0]
                             - pen * main_sequence_masses[0])
    else:
        if mass > main_sequence_masses[-1]:
            tsol = (main_sequence_masses[-1] / mass) * main_sequence_times[-1]
        else:
            index = 0
            while True:
                if mass < main_sequence_masses[index]:
                    pen = ((main_sequence_times[index]
                           - main_sequence_times[index - 1])
                           / (main_sequence_masses[index]
                              - main_sequence_masses[index - 1]))
                    tsol = pen * mass + (main_sequence_times[index]
                                         - pen * main_sequence_masses[index])
                    break

    main_sequence_masses = np.array([0.85, 1.00, 1.25, 1.50, 1.75, 2.00, 3.00])
    # Althaus priv. comm X = 0.752, Y = 0.247
    main_sequence_times = np.array([10.34, 5.756, 2.623, 1.412,
                                    0.905, 0.639, 0.245])

    # TODO: put this in a function as it is the same if as before
    if mass < main_sequence_masses[0]:
        pen = ((main_sequence_times[1] - main_sequence_times[0])
               / (main_sequence_masses[1] - main_sequence_masses[0]))
        tsub = pen * mass + (main_sequence_times[0]
                             - pen * main_sequence_masses[0])
    else:
        if mass > main_sequence_masses[-1]:
            tsub = (main_sequence_masses[-1] / mass) * main_sequence_times[-1]
        else:
            index = 0
            while True:
                if mass < main_sequence_masses[index]:
                    pen = ((main_sequence_times[index]
                           - main_sequence_times[index - 1])
                           / (main_sequence_masses[index]
                              - main_sequence_masses[index - 1]))
                    tsub = pen * mass + (main_sequence_times[index]
                                         - pen * main_sequence_masses[index])
                    break

    return tsub + ((tsol - tsub) / (0.01 - 0.001)) * (metallicity - 0.001)


def get_white_dwarf_masses(progenitor_masses: pd.Series) -> np.ndarray:
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

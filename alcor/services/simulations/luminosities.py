from typing import (Union,
                    Dict,
                    List)

import numpy as np

from alcor.models.star import GalacticDiskType


def get_white_dwarfs(stars: Dict[str, List[Union[float, GalacticDiskType]]],
                     thin_disk_age: float,
                     thick_disk_age: float,
                     halo_age: float,
                     ifmr_parameter: float,
                     chandrasekhar_limit: float = 1.4,
                     max_mass: float = 10.5,
                     solar_metallicity: float = 0.01,
                     subsolar_metallicity: float = 0.001
                     ) -> Dict[str, List[Union[float, GalacticDiskType]]]:
    max_age = max(thin_disk_age, thick_disk_age, halo_age)

    metallicities = []
    cooling_times = []
    masses = []

    for (progenitor_mass,
         galactic_structure_type,
         birth_time) in zip(stars['progenitor_masses'],
                            stars['galactic_structure_types'],
                            stars['birth_times']):
        if progenitor_mass > max_mass:
            continue

        if galactic_structure_type == GalacticDiskType.halo:
            metallicity = subsolar_metallicity
        else:
            metallicity = solar_metallicity

        main_sequence_lifetime = get_main_sequence_life_time(
                mass=progenitor_mass,
                metallicity=metallicity)

        cooling_time = (max_age - birth_time - main_sequence_lifetime)

        if cooling_time < 0.:
            continue

        mass = get_white_dwarf_mass(progenitor_mass) * ifmr_parameter

        if mass > chandrasekhar_limit:
            continue

        metallicities.append(metallicity)
        cooling_times.append(cooling_time)
        masses.append(mass)

    stars['metallicities'] = metallicities
    stars['cooling_times'] = cooling_times
    stars['masses'] = masses

    return stars


# According to model by Leandro & Renedo et al.(2010)
def get_main_sequence_life_time(mass: float,
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


def get_white_dwarf_mass(progenitor_mass: float) -> float:
    if progenitor_mass < 2.7:
        return 0.096 * progenitor_mass + 0.429
    elif 2.7 <= progenitor_mass <= 6.:
        return 0.137 * progenitor_mass + 0.3183
    else:
        return 0.1057 * progenitor_mass + 0.5061

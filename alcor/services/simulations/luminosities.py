from typing import List

import numpy as np

from alcor.models.star import Star


def get_white_dwarfs(stars: List[Star],
                     thin_disk_age: float,
                     ifmr_parameter: float) -> List[Star]:
    white_dwarfs = []

    for star in stars:
        if star.progenitor_mass <= 10.5:
            end_evolution_star = Star()

            end_evolution_star.metallicity = 0.01

            main_sequence_lifetime = get_main_sequence_life_time(
                mass=star.progenitor_mass,
                metallicity=end_evolution_star.metallicity)
            # FIXME: this is only for thin disk stars
            end_evolution_star.cooling_time = (thin_disk_age - star.birth_time
                                               - main_sequence_lifetime)

            if end_evolution_star.cooling_time > 0.:
                end_evolution_star.mass = (get_wd_mass(star.progenitor_mass)
                                           * ifmr_parameter)

            if end_evolution_star.mass <= 1.4:
                white_dwarfs.append(end_evolution_star)

    return white_dwarfs


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


def get_wd_mass(progenitor_mass: float) -> float:
    if progenitor_mass < 2.7:
        return 0.096 * progenitor_mass + 0.429
    elif 2.7 <= progenitor_mass <= 6.:
        return 0.137 * progenitor_mass + 0.3183
    else:
        return 0.1057 * progenitor_mass + 0.5061

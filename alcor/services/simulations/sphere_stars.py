import logging

import numpy as np
import sys

from alcor.models.star import GalacticDiskEnum

logger = logging.getLogger(__name__)


# TODO: take these consts to settings file
def generate_stars_in_sphere(max_stars_count: int = 6E6,
                             time_bins_count: int = 5000,
                             scale_height_factor: float = 1.,
                             thin_disk_scale_height_kpc: float = 0.25,
                             thick_disk_scale_height_kpc: float = 0.9,
                             ttdisk: float = 12.,  # TODO: what is ttdisk?
                             tmdisk: float = 10.,  # TODO: what is tmdisk?
                             tau: float = 2.,  # TODO: what is tau?,
                             stars_birth_init_time_gyr: float = 0.,
                             burst_formation_factor: float = 5.,
                             star_formation_rate_param: float = 25.,
                             sector_radius_kpc: float = 0.05,
                             thin_disk_age_gyr: float = 9.2,
                             burst_age_gyr: float = 0.6,
                             mass_reduction_factor: float = 0.03,
                             thick_disk_stars_fraction: float = 0.8,
                             initial_mass_function_param: float = -2.35
                             ) -> None:
    max_t = tmdisk * np.exp(-tmdisk / tau)  # TODO: what is max_t?

    time_increment = ((thin_disk_age_gyr - stars_birth_init_time_gyr)
                      / time_bins_count)
    sector_area = np.pi * sector_radius_kpc ** 2
    birth_rate = (time_increment * sector_area * 1E6  # TODO: what is 1E6?
                  * mass_reduction_factor
                  * normalization_const(
                        star_formation_rate_param=star_formation_rate_param,
                        thin_disk_age_gyr=thin_disk_age_gyr))
    burst_init_time_gyr = thin_disk_age_gyr - burst_age_gyr

    stars_count = 0

    progenitors_masses = []
    galactic_disk_types = []
    birth_times = []
    z_coordinates = []

    for time_bin in range(time_bins_count):
        total_bin_mass = 0.

        bin_init_time = stars_birth_init_time_gyr + time_bin * time_increment

        if bin_init_time >= burst_init_time_gyr:
            birth_rate *= burst_formation_factor

        while True:
            star_mass = get_mass_from_salpeter_initial_mass_function(
                initial_mass_function_param=initial_mass_function_param)
            stars_count += 1

            if stars_count > max_stars_count:
                logger.error('Number of stars is too high - '
                             'decrease mass reduction factor')
                sys.exit()

            progenitors_masses.append(star_mass)
            total_bin_mass += star_mass

            if np.random.rand() <= thick_disk_stars_fraction:
                galactic_disk_types.append(GalacticDiskEnum.thick)

                while True:
                    # TODO: find out what is going on here
                    ttry = ttdisk * np.random.rand()
                    ft = ttry * np.exp(-ttry / tau)
                    fz = max_t * np.random.rand()
                    if fz <= ft:
                        birth_times.append(ttry)
                        break
            else:
                galactic_disk_types.append(GalacticDiskEnum.thin)
                birth_times.append(stars_birth_init_time_gyr
                                   + time_bin * time_increment
                                   + time_increment * np.random.rand())

            scale_height = (thin_disk_scale_height_kpc
                            if galactic_disk_types[-1] == GalacticDiskEnum.thin
                            else thick_disk_scale_height_kpc)

            # TODO: find out what is going on here
            while True:
                xx = scale_height_factor * scale_height * np.random.rand()
                if xx != 0.:
                    break
            zz = scale_height * np.log(scale_height_factor * scale_height / xx)
            in_param = int(2. * np.random.rand()[0])
            z_coordinates.append(zz * float(1 - 2 * in_param))

            # TODO: find out what is going on here
            if total_bin_mass > birth_rate:
                progenitors_masses[-1] -= total_bin_mass - birth_rate
                break


def normalization_const(star_formation_rate_param: float,
                        thin_disk_age_gyr: float,
                        sigma: float = 51.  # TODO: what is sigma?
                        ) -> float:
    return sigma / (star_formation_rate_param
                    * (1 - np.exp(-thin_disk_age_gyr
                                  / star_formation_rate_param)))


# TODO: implement inverse transform sampling
def get_mass_from_salpeter_initial_mass_function(
        initial_mass_function_param: float,
        min_mass: float = 0.4,
        max_mass: float = 50.
        ) -> float:
    y_max = min_mass ** initial_mass_function_param

    mass_amplitude = max_mass - min_mass
    while True:
        # TODO: implement seeds tracking
        y = y_max * np.random.rand()
        mass = min_mass + mass_amplitude * np.random.rand()[0]
        y_imf = mass ** initial_mass_function_param
        if y <= y_imf:
            return mass

from typing import List, Iterable

import numpy as np

from alcor.models.star import Star, GalacticDiskEnum


# TODO: take all consts to settings file
def generate_stars(*,
                   delta_latitude: float = np.deg2rad(2.64),
                   normalization_cone_height_kpc: float = 0.2,
                   cone_height_kpc: float = 2.,
                   thin_disk_density: float = 0.095 * 1E9,
                   thin_disk_scale_height_kpc: float = 0.25,
                   thick_disk_scale_height_kpc: float = 1.5,
                   thick_disk_scale_length_kpc: float = 3.,
                   solar_galactocentric_distance_kpc: float = 8.5,
                   tmdisk: float = 10.,  # TODO: what is tmdisk?
                   ttdisk: float = 12.,  # TODO: what is ttdisk?
                   tau: float = 2.,  # TODO: what is tau?
                   thin_disk_age: float,
                   mass_reduction_factor: float,
                   thick_disk_stars_fraction: float,
                   cone_height_longitude: float,
                   cone_height_latitude: float,
                   initial_mass_function_param: float = -2.35) -> List[Star]:
    delta_longitude = delta_latitude / np.cos(cone_height_latitude)
    min_longitude = cone_height_longitude - delta_longitude / 2.
    min_latitude = cone_height_latitude - delta_latitude / 2.

    thin_disk_stars_fraction = 1. - thick_disk_stars_fraction

    thin_disk_stars = []
    thick_disk_stars = []

    if thin_disk_stars_fraction > 0.:
        density = thin_disk_density * thin_disk_stars_fraction
        thin_disk_stars = generate(
            cone_height_longitude=cone_height_longitude,
            cone_height_latitude=cone_height_latitude,
            min_longitude=min_longitude,
            min_latitude=min_latitude,
            delta_longitude=delta_longitude,
            cone_height=cone_height_kpc,
            delta_latitude=delta_latitude,
            density=density,
            scale_height=thin_disk_scale_height_kpc,
            normalization_cone_height=normalization_cone_height_kpc,
            mass_reduction_factor=mass_reduction_factor,
            solar_galactocentric_distance=solar_galactocentric_distance_kpc,
            thick_disk_scale_length=thick_disk_scale_length_kpc,
            initial_mass_function_param=initial_mass_function_param)

        for star in thin_disk_stars:
            star.disk_belonging = GalacticDiskEnum.thin
            star.birth_time = thin_disk_age * np.random.rand()[0]

    if thick_disk_stars_fraction > 0.:
        density = thin_disk_density * thick_disk_stars_fraction
        thick_disk_stars = generate(
            cone_height_longitude=cone_height_longitude,
            cone_height_latitude=cone_height_latitude,
            min_longitude=min_longitude,
            min_latitude=min_latitude,
            delta_longitude=delta_longitude,
            cone_height=cone_height_kpc,
            delta_latitude=delta_latitude,
            density=density,
            scale_height=thick_disk_scale_height_kpc,
            normalization_cone_height=normalization_cone_height_kpc,
            mass_reduction_factor=mass_reduction_factor,
            solar_galactocentric_distance=solar_galactocentric_distance_kpc,
            thick_disk_scale_length=thick_disk_scale_length_kpc,
            initial_mass_function_param=initial_mass_function_param)

        for star in thick_disk_stars:
            star.disk_belonging = GalacticDiskEnum.thick
            star.birth_time = get_thick_disk_star_birth_time(tmdisk=tmdisk,
                                                             ttdisk=ttdisk,
                                                             tau=tau)

    return thin_disk_stars + thick_disk_stars


def generate(cone_height_longitude: float,
             cone_height_latitude: float,
             min_longitude: float,
             min_latitude: float,
             delta_longitude: float,
             cone_height: float,
             delta_latitude: float,
             density: float,
             scale_height: float,
             normalization_cone_height: float,
             mass_reduction_factor: float,
             solar_galactocentric_distance: float,
             thick_disk_scale_length: float,
             initial_mass_function_param: float
             ) -> Iterable[Star]:
    normalization_cone_mass = (cone_mass(latitude=cone_height_longitude,
                                         delta_latitude=delta_latitude,
                                         density=density,
                                         scale_height=scale_height,
                                         cone_height=normalization_cone_height)
                               * mass_reduction_factor)
    max_density = get_max_density(
        longitude=cone_height_longitude,
        latitude=cone_height_latitude,
        scale_height=scale_height,
        cone_height=cone_height,
        solar_galactocentric_distance=solar_galactocentric_distance,
        thick_disk_scale_length=thick_disk_scale_length)

    normalization_mass = 0.

    while True:
        star = Star()
        longitude = (min_longitude + delta_longitude * np.random.rand()[0])
        latitude = (min_latitude + delta_latitude * np.random.rand()[0])

        while True:
            distance = cone_height * np.random.rand()[0]
            density = get_density(
                distance=distance,
                longitude=longitude,
                latitude=latitude,
                scale_height=scale_height,
                solar_galactocentric_distance=solar_galactocentric_distance,
                thick_disk_scale_length=thick_disk_scale_length)

            # Monte-Carlo accept/reject method
            random_valid_density = max_density * np.random.rand()[0]
            if random_valid_density <= density:
                break

        # TODO: this get_mass is also in another PR
        star.progenitor_mass = get_mass(initial_mass_function_param)
        # TODO: I am not sure if I need these
        star.r_cylindric_coordinate = opposite_triangle_side(
            solar_galactocentric_distance,
            distance * np.abs(np.cos(latitude)),
            longitude)
        star.th_cylindric_coordinate = (
            np.arcsin(distance * np.abs(np.cos(latitude)))
            * np.sin(longitude) / star.r_cylindric_coordinate)
        star.z_cylindric_coordinate = distance * np.sin(latitude)

        if distance < normalization_cone_height:
            normalization_mass += star.progenitor_mass

        if normalization_mass > normalization_cone_mass:
            break

        yield star


def cone_mass(latitude: float,
              delta_latitude: float,
              density: float,
              scale_height: float,
              cone_height: float,
              max_safe_latitude: float = np.deg2rad(85.)) -> float:
    # Stars distribution is symmetrical with respect to the galactic plane
    latitude = np.abs(latitude)

    # To prevent problems with spherical coordinates near poles
    latitude = min(latitude, max_safe_latitude)

    delta_longitude = delta_latitude / np.cos(latitude)
    min_latitude = latitude - delta_latitude / 2.
    max_latitude = latitude + delta_latitude / 2.

    # Next, we consider 3 cases:
    #     1) both min_latitude and max_latitude are in [0; pi/2]
    #     2) `min_latitude` is negative and `max_latitude` is in [0; pi/2]
    #     3) `min_latitude` is in [0; pi/2] and `max_latitude` > pi/2
    
    # case 2
    if min_latitude < 0.:
        mass = (density * delta_longitude
                * (iota_integral(np.abs(min_latitude),
                                 scale_height,
                                 cone_height)
                   - iota_integral(max_latitude,
                                   scale_height,
                                   cone_height)))
    # case 3
    elif max_latitude > np.pi / 2.:
        mass = (density * delta_longitude
                * (lambda_integral(np.pi - max_latitude,
                                   scale_height,
                                   cone_height)
                   + lambda_integral(min_latitude,
                                     scale_height,
                                     cone_height)))
    # case 1
    else:
        mass = (density * delta_longitude * scale_height
                * (kappa_integral(min_latitude, scale_height, cone_height)
                   - kappa_integral(max_latitude, scale_height, cone_height)))
    return mass


def kappa_integral(latitude: float,
                   scale_height: float,
                   height: float) -> float:
    # Next var is used just to simplify the following expression
    gamma = scale_height / np.sin(latitude)
    return gamma * (gamma - np.exp(-height / gamma) * (gamma + height))


def iota_integral(latitude: float,
                  scale_height: float,
                  height: float) -> float:
    return (height ** 2 / 2.
            - scale_height * kappa_integral(latitude, scale_height, height))


def lambda_integral(latitude: float,
                    scale_height: float,
                    height: float) -> float:
    return (-scale_height * (kappa_integral(np.pi / 2., scale_height, height)
                             - kappa_integral(latitude, scale_height, height)))


def get_max_density(*,
                    longitude: float,
                    latitude: float,
                    scale_height: float,
                    cone_height: float,
                    distance_bins_count: int = 1000,
                    monte_carlo_shift_factor: float = 1.1,
                    solar_galactocentric_distance: float,
                    thick_disk_scale_length: float
                    ) -> float:
    distance_bin_size = cone_height / distance_bins_count

    max_density = 0.

    for distance_bin in range(distance_bins_count):
        distance = distance_bin_size * (distance_bin + 1)
        density = get_density(
            distance=distance,
            longitude=longitude,
            latitude=latitude,
            scale_height=scale_height,
            solar_galactocentric_distance=solar_galactocentric_distance,
            thick_disk_scale_length=thick_disk_scale_length)
        max_density = np.maximum(max_density, density)

    max_density *= monte_carlo_shift_factor
    return max_density


def get_density(distance: float,
                longitude: float,
                latitude: float,
                scale_height: float,
                solar_galactocentric_distance: float,
                thick_disk_scale_length: float) -> float:
    pole_projection = distance * np.abs(np.sin(latitude))
    plane_projection = distance * np.abs(np.cos(latitude))

    galactocentric_distance = opposite_triangle_side(
        plane_projection,
        solar_galactocentric_distance,
        longitude)

    density = (distance ** 2 * np.exp(-np.abs(pole_projection) / scale_height)
               * np.exp(-np.abs(galactocentric_distance)
                        / thick_disk_scale_length))
    return density


def opposite_triangle_side(adjacent_1: float,
                           adjacent_2: float,
                           enclosed_angle: float) -> float:
    return np.sqrt(adjacent_1 ** 2 + adjacent_2 ** 2
                   - 2. * adjacent_1 * adjacent_2 * np.cos(enclosed_angle))


def get_mass(initial_mass_function_param: float,
             min_mass: float = 0.4,
             max_mass: float = 50.) -> float:
    mass_range = max_mass - min_mass
    y_max = min_mass ** initial_mass_function_param

    while True:
        y = y_max * np.random.rand()[0]
        mass = min_mass + mass_range * np.random.rand()[0]
        y_imf = mass ** initial_mass_function_param

        if y <= y_imf:
            break

    return mass


def get_thick_disk_star_birth_time(tmdisk: float,
                                   ttdisk: float,
                                   tau: float) -> float:
    max_t = tmdisk * np.exp(-tmdisk / tau)
    while True:
        ttry = ttdisk * np.random.rand()[0]
        ft = ttry * np.exp(-ttry / tau)
        fz = max_t * np.random.rand()[0]
        if fz <= ft:
            return ttry

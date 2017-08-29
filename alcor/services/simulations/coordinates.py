from typing import List

import numpy as np

from alcor.models.star import Star


def calculate_coordinates(stars: List[Star],
                          solar_galactocentric_distance: float,
                          deltag: float = 0.478,  # TODO: what is deltag?
                          theta: float = 2.147,  # TODO: what is theta?
                          alphag: float = 3.35  # TODO: what is alphag
                          ) -> None:
    for star in stars:
        # TODO: give shorter name
        stellar_galactocentric_distance_plane_projection = (
            opposite_triangle_side(solar_galactocentric_distance,
                                   star.r_cylindric_coordinate,
                                   star.th_cylindric_coordinate))
        star.galactic_distance = np.sqrt(
            stellar_galactocentric_distance_plane_projection ** 2
            + star.z_coordinate ** 2)

        # TODO: implement function
        star.galactic_longitude = np.arccos(
            (solar_galactocentric_distance ** 2
             + stellar_galactocentric_distance_plane_projection ** 2
             - star.r_cylindric_coordinate ** 2)
            / (2. * stellar_galactocentric_distance_plane_projection
               * solar_galactocentric_distance))

        if (star.r_cylindric_coordinate * np.cos(star.th_cylindric_coordinate)
                > solar_galactocentric_distance):
            star.galactic_longitude = np.pi - star.galactic_longitude
        elif np.sin(star.th_cylindric_coordinate) < 0.:
            star.galactic_longitude += 2. * np.pi

        if star.galactic_longitude > 2. * np.pi:
            star.galactic_longitude -= 2. * np.pi

        # TODO: or use arctan2
        star.galactic_latitude = np.arctan(
            np.abs(star.z_coordinate
                   / stellar_galactocentric_distance_plane_projection))

        if star.z_coordinate < 0.:
            star.galactic_latitude = -star.galactic_latitude

        sin_longitude = np.sin(star.galactic_longitude)
        cos_longitude = np.cos(star.galactic_longitude)
        sin_latitude = np.sin(star.galactic_latitude)
        cos_latitude = np.cos(star.galactic_latitude)

        # TODO: what is this?
        zkri = 1. / (4.74E3 * star.galactic_distance)

        star.proper_motion_component_l = (
            (-zkri * (sin_longitude / cos_latitude) * star.u_velocity)
            + (zkri * (cos_longitude / cos_latitude) * star.v_velocity))
        star.proper_motion_component_b = (
            (-zkri * cos_longitude * sin_latitude * star.u_velocity)
            + (-zkri * sin_latitude * sin_longitude * star.v_velocity)
            + (zkri * cos_latitude * star.w_velocity))
        # TODO: rename as radial velocity
        star.proper_motion_component_vr = (
            (cos_latitude * cos_longitude * star.u_velocity)
            + (cos_latitude * sin_latitude * star.v_velocity)
            + (sin_latitude * star.w_velocity))
        star.proper_motion = np.sqrt(star.proper_motion_component_l ** 2
                                     + star.proper_motion_component_b ** 2)

        star.declination = (np.arcsin(np.sin(deltag) * sin_latitude
                            + np.cos(deltag) * cos_latitude
                              * np.cos(theta - star.galactic_longitude)))

        # TODO: what is xs and xc?
        xs = ((cos_latitude * np.sin(theta - star.galactic_longitude))
              / np.cos(star.declination))
        xc = ((np.cos(deltag) * sin_latitude - np.sin(deltag) * cos_latitude
               * np.cos(theta - star.galactic_longitude))
              / np.cos(star.declination))

        if xs >= 0.:
            if xc >= 0.:
                star.right_ascension = np.arcsin(xs) + alphag
            else:
                star.right_ascension = np.arccos(xc) + alphag
        else:
            if xc < 0.:
                star.right_ascension = np.pi - np.arcsin(xs) + alphag
            else:
                star.right_ascension = 2. * np.pi + np.arcsin(xs) + alphag

        if star.right_ascension > 2. * np.pi:
            star.right_ascension -= 2. * np.pi


def opposite_triangle_side(adjacent_1: float,
                           adjacent_2: float,
                           enclosed_angle: float) -> float:
    return np.sqrt(adjacent_1 ** 2 + adjacent_2 ** 2
                   - 2. * adjacent_1 * adjacent_2 * np.cos(enclosed_angle))
from typing import Union

import numpy as np
import pandas as pd

PC_PER_KPC = 1e3


def set_coordinates(stars: pd.DataFrame,
                    *,
                    solar_galactocentric_distance: float) -> None:
    """
    More info on conversions at:
    https://physics.stackexchange.com/questions/88663/converting-between-galactic-and-ecliptic-coordinates
    :param stars: stars without coordinates
    :param solar_galactocentric_distance: distance from Sun to Galactic Center
    :return: stars with assigned coordinates and proper motions
    """
    distance_plane_projections = (opposite_triangle_side(
            solar_galactocentric_distance,
            stars['r_cylindrical'].values,
            stars['theta_cylindrical'].values))
    stars['distance'] = np.sqrt(distance_plane_projections ** 2
                                + stars['z_coordinate'] ** 2)
    stars['galactic_longitude'] = get_galactic_longitudes(
            solar_galactocentric_distance=solar_galactocentric_distance,
            r_cylindrical=stars['r_cylindrical'].values,
            thetas_cylindrical=stars['theta_cylindrical'].values,
            distance_plane_projections=distance_plane_projections)
    stars['galactic_latitude'] = get_galactic_latitudes(
            z_coordinates=stars['z_coordinate'].values,
            distance_plane_projections=distance_plane_projections)


def get_galactic_latitudes(*,
                           z_coordinates: np.ndarray,
                           distance_plane_projections: np.ndarray
                           ) -> np.ndarray:
    # TODO: or use arctan2?
    result = np.arctan(np.abs(z_coordinates / distance_plane_projections))
    return np.copysign(result, z_coordinates)


def get_galactic_longitudes(*,
                            solar_galactocentric_distance: float,
                            r_cylindrical: np.ndarray,
                            thetas_cylindrical: np.ndarray,
                            distance_plane_projections: np.ndarray
                            ) -> np.ndarray:
    result = triangle_angle(adjacent=np.array([solar_galactocentric_distance]),
                            other_adjacent=distance_plane_projections,
                            opposite=r_cylindrical)
    # Unfolding from 0-180 to 0-360
    mask = thetas_cylindrical > np.pi
    result[mask] = 2. * np.pi - result[mask]

    return result


def opposite_triangle_side(adjacent: Union[float, np.ndarray],
                           other_adjacent: Union[float, np.ndarray],
                           enclosed_angle: Union[float, np.ndarray]
                           ) -> Union[float, np.ndarray]:
    return np.sqrt(adjacent ** 2 + other_adjacent ** 2
                   - 2. * adjacent * other_adjacent * np.cos(enclosed_angle))


def triangle_angle(adjacent: np.ndarray,
                   other_adjacent: np.ndarray,
                   opposite: np.ndarray) -> np.ndarray:
    arccos_argument = ((adjacent ** 2 + other_adjacent ** 2 - opposite ** 2)
                       / (2. * adjacent * other_adjacent))

    positive_erroneous_values_mask = arccos_argument > 1.
    arccos_argument[positive_erroneous_values_mask] = 1.

    negative_erroneous_values_mask = arccos_argument < -1.
    arccos_argument[negative_erroneous_values_mask] = -1.

    return np.arccos(arccos_argument)

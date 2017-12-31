import numpy as np

from alcor.services.simulations.coordinates import (triangle_angle,
                                                    opposite_triangle_side,
                                                    get_galactic_longitudes,
                                                    get_galactic_latitudes,
                                                    right_ascensions)


def test_triangle_angle(adjacent: float,
                        other_adjacent: float,
                        opposite: float) -> None:
    angle = triangle_angle(np.array([adjacent]),
                           np.array([other_adjacent]),
                           np.array([opposite]))

    assert isinstance(angle, np.ndarray)
    assert angle >= 0.


def test_opposite_triangle_side(adjacent: float,
                                other_adjacent: float,
                                enclosed_angle: float) -> None:
    side = opposite_triangle_side(adjacent=np.array([adjacent]),
                                  other_adjacent=np.array([other_adjacent]),
                                  enclosed_angle=np.array([enclosed_angle]))

    assert isinstance(side, np.ndarray)
    assert np.asscalar(side) > 0.


def test_get_galactic_longitudes(solar_galactocentric_distance: float,
                                 single_r_cylindrical: float,
                                 theta_cylindrical: float,
                                 distance_plane_projections: float) -> None:
    longitude = get_galactic_longitudes(
            solar_galactocentric_distance=solar_galactocentric_distance,
            r_cylindrical=np.array([single_r_cylindrical]),
            thetas_cylindrical=np.array([theta_cylindrical]),
            distance_plane_projections=np.array([distance_plane_projections]))
    assert isinstance(longitude, np.ndarray)
    assert np.pi * 2. > np.asscalar(longitude) >= 0.


def test_get_galactic_latitudes(z_coordinates: float,
                                distance_plane_projections: float) -> None:
    latitude = get_galactic_latitudes(
            z_coordinates=np.array([z_coordinates]),
            distance_plane_projections=np.array([distance_plane_projections]))
    assert isinstance(latitude, np.ndarray)
    assert np.pi * 2. > np.asscalar(latitude) >= 0.


def test_right_ascensions(cos_latitude: float,
                          sin_latitude: float,
                          theta: float,
                          galactic_longitude: float,
                          declination: float,
                          ngp_declination: float,
                          ngp_right_ascension: float) -> None:
    angle = right_ascensions(
            cos_latitude=np.array([cos_latitude]),
            sin_latitude=np.array([sin_latitude]),
            theta=theta,
            galactic_longitudes=np.array([galactic_longitude]),
            declinations=np.array([declination]),
            ngp_declination=ngp_declination,
            ngp_right_ascension=ngp_right_ascension)

    assert isinstance(angle, np.ndarray)
    assert np.pi * 2. > np.asscalar(angle) >= 0.

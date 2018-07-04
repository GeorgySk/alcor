import numpy as np
import pandas as pd


def assign_equatorial_coordinates(stars: pd.DataFrame,
                                  *,
                                  ngp_declination: float = 3.35,
                                  ngp_right_ascension: float = 3.35,
                                  theta: float = 2.147,
                                  sin_latitude: np.ndarray,
                                  cos_latitude: np.ndarray):
    """
    :param stars:
    :param ngp_declination: declination of the North Galactic Pole
    :param ngp_right_ascension: right ascension of the North Galactic Pole
    :param theta: see link above, same as BK = 122.9º there
    :param sin_latitude:
    :param cos_latitude:
    :return:
    """
    stars['declination'] = (np.arcsin(
            np.sin(ngp_declination) * sin_latitude
            + np.cos(ngp_declination) * cos_latitude
            * np.cos(theta - stars['galactic_longitude'])))

    stars['right_ascension'] = right_ascensions(
            cos_latitude=cos_latitude,
            sin_latitude=sin_latitude,
            theta=theta,
            galactic_longitudes=stars['galactic_longitude'],
            declinations=stars['declination'],
            ngp_declination=ngp_declination,
            ngp_right_ascension=ngp_right_ascension)


def right_ascensions(*,
                     cos_latitude: np.ndarray,
                     sin_latitude: np.ndarray,
                     theta: float,
                     galactic_longitudes: np.ndarray,
                     declinations: np.ndarray,
                     ngp_declination: float,
                     ngp_right_ascension: float) -> np.ndarray:
    result = np.empty(galactic_longitudes.size)

    # TODO: give better names
    # These variables are for conversion from galactic to equatorial
    # coordinates. More info at the link above
    xs = (cos_latitude * np.sin(theta - galactic_longitudes)
          / np.cos(declinations))
    xc = ((np.cos(ngp_declination) * sin_latitude
           - np.sin(ngp_declination) * cos_latitude
             * np.cos(theta - galactic_longitudes)) / np.cos(declinations))

    # TODO: find out what is going on here
    xs_ge_zero_xc_ge_zero_mask = ((xs >= 0.) & (xc >= 0.))
    xs_ge_zero_xc_lt_zero_mask = ((xs >= 0.) & (xc < 0.))
    xs_lt_zero_xc_ge_zero_mask = ((xs < 0.) & (xc >= 0.))
    xs_lt_zero_xc_lt_zero_mask = ((xs < 0.) & (xc < 0.))

    result[xs_ge_zero_xc_ge_zero_mask] = (
        np.arcsin(xs[xs_ge_zero_xc_ge_zero_mask]) + ngp_right_ascension)
    result[xs_ge_zero_xc_lt_zero_mask] = (
        np.arccos(xc[xs_ge_zero_xc_lt_zero_mask]) + ngp_right_ascension)
    result[xs_lt_zero_xc_ge_zero_mask] = (
        2. * np.pi + np.arcsin(xs[xs_lt_zero_xc_ge_zero_mask])
        + ngp_right_ascension)
    result[xs_lt_zero_xc_lt_zero_mask] = (
        np.pi - np.arcsin(xs[xs_lt_zero_xc_lt_zero_mask])
        + ngp_right_ascension)

    mask = result > 2. * np.pi
    result[mask] -= 2. * np.pi

    return result

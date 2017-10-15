import numpy as np
import pandas as pd


def by_parallax(stars: pd.DataFrame,
                *,
                min_parallax: float) -> pd.DataFrame:
    distances_in_pc = stars['distance'] * 1e3
    parallaxes = 1 / distances_in_pc
    return stars[parallaxes > min_parallax]


def by_declination(stars: pd.DataFrame,
                   *,
                   min_declination: float) -> pd.DataFrame:
    return stars[stars['declination'] > min_declination]


def by_velocity(stars: pd.DataFrame,
                *,
                max_velocity: float) -> pd.DataFrame:
    return stars[np.power(stars['u_velocity'], 2)
                 + np.power(stars['v_velocity'], 2)
                 + np.power(stars['w_velocity'], 2)
                 < max_velocity ** 2]


def by_proper_motion(stars: pd.DataFrame,
                     *,
                     min_proper_motion: float) -> pd.DataFrame:
    return stars[stars['proper_motion'] > min_proper_motion]


# TODO: find out what is going on here
def by_reduced_proper_motion(stars: pd.DataFrame) -> pd.DataFrame:
    # Transformation from UBVRI to ugriz. More info at:
    # Jordi, Grebel & Ammon, 2006, A&A, 460; equations 1-8 and Table 3
    g_ugriz_abs_magnitudes = (stars['v_abs_magnitude'] - 0.124
                              + 0.63 * (stars['b_abs_magnitude']
                                        - stars['v_abs_magnitude']))
    z_ugriz_abs_magnitudes = (g_ugriz_abs_magnitudes
                              - 1.646 * (stars['v_abs_magnitude']
                                         - stars['r_abs_magnitude'])
                              - 1.584 * (stars['r_abs_magnitude']
                                         - stars['i_abs_magnitude'])
                              + 0.525)
    g_apparent_magnitudes = apparent_magnitude(g_ugriz_abs_magnitudes,
                                               distance_kpc=stars['distance'])
    z_apparent_magnitudes = apparent_magnitude(z_ugriz_abs_magnitudes,
                                               distance_kpc=stars['distance'])
    # TODO: find out the meaning and check if the last 5 is correct
    hrms = g_apparent_magnitudes + 5. * np.log10(stars['proper_motion']) + 5.
    stars = stars[(g_apparent_magnitudes - z_apparent_magnitudes > -0.33)
                  | (hrms > 14.)]

    return stars[hrms > 15.17 + 3.559 * (g_apparent_magnitudes
                                         - z_apparent_magnitudes)]


def by_apparent_magnitude(stars: pd.DataFrame,
                          *,
                          max_v_apparent_magnitude: float) -> pd.DataFrame:
    v_apparent_magnitudes = apparent_magnitude(stars['v_abs_magnitude'],
                                               distance_kpc=stars['distance'])
    return stars[v_apparent_magnitudes <= max_v_apparent_magnitude]


def apparent_magnitude(abs_magnitude: pd.Series,
                       distance_kpc: pd.Series) -> pd.Series:
    # More info at (2nd formula, + 3.0 because the distance is in kpc):
    # https://en.wikipedia.org/wiki/Absolute_magnitude#Apparent_magnitude
    return abs_magnitude - 5. + 5. * (np.log10(distance_kpc) + 3.)

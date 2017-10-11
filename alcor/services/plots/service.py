import uuid
from collections import Counter
from functools import partial
from typing import Set

import numpy as np
import pandas as pd
from sqlalchemy.orm.session import Session
from sqlalchemy.sql import text

from alcor.models import eliminations
from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds,
               heatmaps,
               toomre_diagram,
               ugriz_diagrams)

ASTRONOMICAL_UNIT = 4.74


def draw(*,
         group_id: uuid.UUID,
         filtration_method: str,
         nullify_radial_velocity: bool,
         with_luminosity_function: bool,
         with_velocities_vs_magnitude: bool,
         with_velocity_clouds: bool,
         lepine_criterion: bool,
         heatmaps_axes: str,
         with_toomre_diagram: bool,
         with_ugriz_diagrams: bool,
         desired_stars_count: int,
         session: Session) -> None:
    fields_to_fetch = get_fields_to_fetch(
            filtration_method=filtration_method,
            nullify_radial_velocity=nullify_radial_velocity,
            lepine_criterion=lepine_criterion,
            with_luminosity_function=with_luminosity_function,
            with_velocities_vs_magnitude=with_velocities_vs_magnitude,
            with_velocity_clouds=with_velocity_clouds,
            heatmaps_axes=heatmaps_axes,
            with_toomre_diagram=with_toomre_diagram,
            with_ugriz_diagrams=with_ugriz_diagrams)

    if not fields_to_fetch:
        return

    fields_str = ", ".join(fields_to_fetch)

    # TODO: find out if it's possible to do without text
    sql = text(f"SELECT {fields_str} FROM stars "
               f"WHERE group_id='{group_id}'")

    if desired_stars_count:
        sql = text(str(sql)
                   + f" ORDER BY RANDOM() LIMIT {desired_stars_count}")

    stars = pd.read_sql_query(sql=sql,
                              con=session.get_bind(),
                              index_col='id')

    stars = filter_stars(stars,
                         method=filtration_method,
                         group_id=group_id,
                         session=session)

    if nullify_radial_velocity:
        stars = set_radial_velocity_to_zero(stars)

    if with_luminosity_function:
        luminosity_function.plot(stars=stars)

    if with_velocities_vs_magnitude:
        if lepine_criterion:
            velocities_vs_magnitude.plot_lepine_case(stars=stars)
        else:
            velocities_vs_magnitude.plot(stars=stars)

    if with_velocity_clouds:
        if lepine_criterion:
            velocity_clouds.plot_lepine_case(stars=stars)
        else:
            velocity_clouds.plot(stars=stars)

    if heatmaps_axes:
        heatmaps.plot(stars=stars,
                      axes=heatmaps_axes)

    if with_toomre_diagram:
        toomre_diagram.plot(stars=stars)

    if with_ugriz_diagrams:
        ugriz_diagrams.plot(stars=stars)


def get_fields_to_fetch(filtration_method: str,
                        nullify_radial_velocity: bool,
                        lepine_criterion: bool,
                        with_luminosity_function: bool,
                        with_velocities_vs_magnitude: bool,
                        with_velocity_clouds: bool,
                        heatmaps_axes: str,
                        with_toomre_diagram: bool,
                        with_ugriz_diagrams: bool) -> Set:
    fields_to_fetch = set()

    if filtration_method != 'raw':
        fields_to_fetch.update(['distance',
                                'declination',
                                'u_velocity',
                                'v_velocity',
                                'w_velocity'])
    if filtration_method == 'restricted':
        fields_to_fetch.update(['b_abs_magnitude',
                                'v_abs_magnitude',
                                'r_abs_magnitude',
                                'i_abs_magnitude',
                                'proper_motion'])

    if nullify_radial_velocity:
        fields_to_fetch.update(['galactic_longitude',
                                'galactic_latitude',
                                'proper_motion_component_l',
                                'proper_motion_component_b',
                                'distance'])

    if lepine_criterion:
        fields_to_fetch.update(['right_ascension',
                                'declination',
                                'distance'])

    if with_luminosity_function:
        fields_to_fetch.update(['luminosity'])

    if with_velocities_vs_magnitude:
        fields_to_fetch.update(['luminosity',
                                'u_velocity',
                                'v_velocity',
                                'w_velocity'])

    if heatmaps_axes or with_velocity_clouds:
        fields_to_fetch.update(['u_velocity',
                                'v_velocity',
                                'w_velocity'])

    if with_toomre_diagram:
        fields_to_fetch.update(['u_velocity',
                                'v_velocity',
                                'w_velocity',
                                'spectral_type'])

    if with_ugriz_diagrams:
        fields_to_fetch.update(['b_abs_magnitude',
                                'v_abs_magnitude',
                                'r_abs_magnitude',
                                'i_abs_magnitude',
                                'spectral_type'])

    if fields_to_fetch:
        fields_to_fetch.update(['id'])

    return fields_to_fetch


def apparent_magnitude(abs_magnitude: pd.Series,
                       distance_kpc: pd.Series
                       ) -> pd.Series:
    # More info at (2nd formula, + 3.0 because the distance is in kpc):
    # https://en.wikipedia.org/wiki/Absolute_magnitude#Apparent_magnitude
    return abs_magnitude - 5. + 5. * (np.log10(distance_kpc) + 3.)


def set_radial_velocity_to_zero(stars: pd.DataFrame) -> pd.DataFrame:
    distances_in_pc = stars['distance'] * 1e3

    a1 = (-ASTRONOMICAL_UNIT * np.cos(stars['galactic_latitude'])
          * np.sin(stars['galactic_longitude']))
    b1 = (-ASTRONOMICAL_UNIT * np.sin(stars['galactic_latitude'])
          * np.cos(stars['galactic_longitude']))
    stars['u_velocity'] = ((a1 * stars['proper_motion_component_l']
                            + b1 * stars['proper_motion_component_b'])
                              * distances_in_pc)

    a2 = (ASTRONOMICAL_UNIT * np.cos(stars['galactic_latitude'])
          * np.cos(stars['galactic_longitude']))
    b2 = (-ASTRONOMICAL_UNIT * np.sin(stars['galactic_latitude'])
          * np.sin(stars['galactic_longitude']))
    stars['v_velocity'] = ((a2 * stars['proper_motion_component_l']
                            + b2 * stars['proper_motion_component_b'])
                              * distances_in_pc)

    b3 = ASTRONOMICAL_UNIT * np.cos(stars['galactic_latitude'])
    stars['w_velocity'] = (b3 * stars['proper_motion_component_b']
                           * distances_in_pc)
    return stars


def filter_stars(stars: pd.DataFrame,
                 *,
                 method: str,
                 group_id: uuid.UUID,
                 session: Session,
                 min_parallax: float = 0.025,
                 min_declination: float = 0.,
                 max_velocity: float = 500.,
                 min_proper_motion: float = 0.04,
                 max_v_apparent_magnitude: float = 19.) -> pd.DataFrame:
    eliminations_counter = Counter()
    eliminations_counter['raw'] = stars.shape[0]

    filtration_functions = {}
    # TODO: fix geometry of a simulated region so that we don't need to use
    # the 'full' filtration method
    if method != 'raw':
        by_parallax = partial(filter_by_parallax,
                              min_parallax=min_parallax)
        by_declination = partial(filter_by_declination,
                                 min_declination=min_declination)
        by_velocity = partial(filter_by_velocity,
                              max_velocity=max_velocity)

        filtration_functions.update(by_parallax=by_parallax,
                                    by_declination=by_declination,
                                    by_velocity=by_velocity)

    if method == 'restricted':
        by_proper_motion = partial(filter_by_proper_motion,
                                   min_proper_motion=min_proper_motion)
        by_apparent_magnitude = partial(
                filter_by_apparent_magnitude,
                max_v_apparent_magnitude=max_v_apparent_magnitude)

        filtration_functions.update(
                by_proper_motion=by_proper_motion,
                by_reduced_proper_motion=filter_by_reduced_proper_motion,
                by_apparent_magnitude=by_apparent_magnitude)

    for criterion, filtration_function in filtration_functions.items():
        stars_count = stars.shape[0]
        stars = filtration_function(stars)
        eliminations_counter[criterion] = stars_count - stars.shape[0]

    counter = eliminations.StarsCounter(group_id=group_id,
                                        **eliminations_counter)
    session.add(counter)
    session.commit()

    return stars


def filter_by_parallax(stars: pd.DataFrame,
                       *,
                       min_parallax: float) -> pd.DataFrame:
    distances_in_pc = stars['distance'] * 1e3
    parallaxes = 1 / distances_in_pc
    return stars[parallaxes > min_parallax]


def filter_by_declination(stars: pd.DataFrame,
                          *,
                          min_declination: float) -> pd.DataFrame:
    return stars[stars['declination'] > min_declination]


def filter_by_velocity(stars: pd.DataFrame,
                       *,
                       max_velocity: float) -> pd.DataFrame:
    return stars[np.power(stars['u_velocity'], 2)
                 + np.power(stars['v_velocity'], 2)
                 + np.power(stars['w_velocity'], 2)
                 < max_velocity ** 2]


def filter_by_proper_motion(stars: pd.DataFrame,
                            *,
                            min_proper_motion: float) -> pd.DataFrame:
    return stars[stars['proper_motion'] > min_proper_motion]


# TODO: find out what is going on here
def filter_by_reduced_proper_motion(stars: pd.DataFrame) -> pd.DataFrame:
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


def filter_by_apparent_magnitude(stars: pd.DataFrame,
                                 *,
                                 max_v_apparent_magnitude: float
                                 ) -> pd.DataFrame:
    v_apparent_magnitudes = apparent_magnitude(stars['v_abs_magnitude'],
                                               distance_kpc=stars['distance'])
    return stars[v_apparent_magnitudes <= max_v_apparent_magnitude]

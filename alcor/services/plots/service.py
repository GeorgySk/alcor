import uuid
from collections import Counter
from typing import Set

import numpy as np
import pandas as pd
from sqlalchemy.engine import Engine
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


def draw(group_id: uuid.UUID,
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
         session: Session,
         engine: Engine) -> None:
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
    if desired_stars_count:
        sql = text(f"SELECT {fields_str} FROM stars " 
                   f"WHERE group_id='{group_id}' "
                   f"ORDER BY RANDOM() LIMIT {desired_stars_count}")
    else:
        sql = text(f"SELECT {fields_str} FROM stars " 
                   f"WHERE group_id='{group_id}'")

    stars = pd.read_sql_query(sql=sql,
                              con=engine,
                              index_col='id')

    eliminations_counter = Counter()

    stars_count = stars.shape[0]
    eliminations_counter['raw'] = stars_count

    if filtration_method != 'raw':
        distances_in_pc = stars['distance'] * 1e3
        parallaxes = 1 / distances_in_pc
        min_parallax = 0.025
        stars = stars[parallaxes > min_parallax]
        eliminations_counter['by_parallax'] = stars_count - stars.shape[0]

        stars_count = stars.shape[0]
        min_declination = 0.
        stars = stars[stars['declination'] > min_declination]
        eliminations_counter['by_declination'] = stars_count - stars.shape[0]

        stars_count = stars.shape[0]
        max_velocity = 500.
        stars = stars[np.power(stars['u_velocity'], 2)
                      + np.power(stars['v_velocity'], 2)
                      + np.power(stars['w_velocity'], 2) < max_velocity ** 2]
        eliminations_counter['by_velocity'] = stars_count - stars.shape[0]

    if filtration_method == 'restricted':
        stars_count = stars.shape[0]
        min_proper_motion = 0.04
        stars = stars[stars['proper_motion'] > min_proper_motion]
        eliminations_counter['by_proper_motion'] = stars_count - stars.shape[0]

        stars_count = stars.shape[0]
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
        g_apparent_magnitudes = apparent_magnitude(
                g_ugriz_abs_magnitudes,
                distance_kpc=stars['distance'])
        z_apparent_magnitudes = apparent_magnitude(
                z_ugriz_abs_magnitudes,
                distance_kpc=stars['distance'])
        # TODO: find out the meaning and check if the last 5 is correct
        hrms = g_apparent_magnitudes + (5. * np.log10(stars['proper_motion'])
                                        + 5.)
        stars = stars[(g_apparent_magnitudes - z_apparent_magnitudes > -0.33)
                      | (hrms > 14.)]
        stars = stars[hrms > 3.559 * (g_apparent_magnitudes
                                      - z_apparent_magnitudes) + 15.17]
        eliminations_counter['by_reduced_proper_motion'] = (stars_count -
                                                            stars.shape[0])

        stars_count = stars.shape[0]
        v_apparent_magnitudes = apparent_magnitude(
                stars['v_abs_magnitude'],
                distance_kpc=stars['distance'])
        stars = stars[v_apparent_magnitudes <= 19.]
        eliminations_counter['by_apparent_magnitude'] = (stars_count -
                                                         stars.shape[0])

    counter = eliminations.StarsCounter(group_id=group_id,
                                        **eliminations_counter)
    session.add(counter)
    session.commit()

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
    stars_copy = stars.copy()

    distances_in_pc = stars['distance'] * 1e3

    a1 = (-ASTRONOMICAL_UNIT * np.cos(stars['galactic_latitude'])
          * np.sin(stars['galactic_longitude']))
    b1 = (-ASTRONOMICAL_UNIT * np.sin(stars['galactic_latitude'])
          * np.cos(stars['galactic_longitude']))
    stars_copy['u_velocity'] = ((a1 * stars['proper_motion_component_l']
                                 + b1 * stars['proper_motion_component_b'])
                                * distances_in_pc)

    a2 = (ASTRONOMICAL_UNIT * np.cos(stars['galactic_latitude'])
          * np.cos(stars['galactic_longitude']))
    b2 = (-ASTRONOMICAL_UNIT * np.sin(stars['galactic_latitude'])
          * np.sin(stars['galactic_longitude']))
    stars_copy['v_velocity'] = ((a2 * stars['proper_motion_component_l']
                                 + b2 * stars['proper_motion_component_b'])
                                * distances_in_pc)

    b3 = ASTRONOMICAL_UNIT * np.cos(stars['galactic_latitude'])
    stars_copy['w_velocity'] = (b3 * stars['proper_motion_component_b']
                                * distances_in_pc)
    return stars_copy

import uuid
from collections import Counter
from functools import partial
from typing import List

import numpy as np
import pandas as pd
from sqlalchemy import func
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlalchemy.orm.session import Session

from alcor.models import eliminations
from alcor.models.star import Star
from . import filters
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
    entities = star_query_entities(
            filtration_method=filtration_method,
            nullify_radial_velocity=nullify_radial_velocity,
            lepine_criterion=lepine_criterion,
            with_luminosity_function=with_luminosity_function,
            with_velocities_vs_magnitude=with_velocities_vs_magnitude,
            with_velocity_clouds=with_velocity_clouds,
            heatmaps_axes=heatmaps_axes,
            with_toomre_diagram=with_toomre_diagram,
            with_ugriz_diagrams=with_ugriz_diagrams)

    if not entities:
        return

    query = (session.query(Star)
             .filter(Star.group_id == group_id))

    if desired_stars_count:
        query = (query.order_by(func.random())
                 .limit(desired_stars_count))

    query = query.with_entities(*entities)

    statement = query.statement

    stars = pd.read_sql_query(sql=statement,
                              con=session.get_bind(),
                              index_col='id')
    stars = filter_stars(stars,
                         method=filtration_method,
                         group_id=group_id,
                         session=session)

    if nullify_radial_velocity:
        set_radial_velocity_to_zero(stars)

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


def set_radial_velocity_to_zero(stars: pd.DataFrame) -> None:
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
        by_parallax = partial(filters.filter_by_parallax,
                              min_parallax=min_parallax)
        by_declination = partial(filters.filter_by_declination,
                                 min_declination=min_declination)
        by_velocity = partial(filters.filter_by_velocity,
                              max_velocity=max_velocity)

        filtration_functions.update(by_parallax=by_parallax,
                                    by_declination=by_declination,
                                    by_velocity=by_velocity)

    if method == 'restricted':
        by_proper_motion = partial(filters.filter_by_proper_motion,
                                   min_proper_motion=min_proper_motion)
        by_apparent_magnitude = partial(
                filters.filter_by_apparent_magnitude,
                max_v_apparent_magnitude=max_v_apparent_magnitude)

        filtration_functions.update(
                by_proper_motion=by_proper_motion,
                by_reduced_proper_motion=(
                    filters.filter_by_reduced_proper_motion),
                by_apparent_magnitude=by_apparent_magnitude)

    for criterion, filtration_function in filtration_functions.items():
        stars_count_before_filtration = stars.shape[0]
        stars = filtration_function(stars)
        eliminations_counter[criterion] = (stars_count_before_filtration
                                           - stars.shape[0])

    counter = eliminations.StarsCounter(group_id=group_id,
                                        **eliminations_counter)
    session.add(counter)
    session.commit()

    return stars


def star_query_entities(*,
                        filtration_method: str,
                        nullify_radial_velocity: bool,
                        lepine_criterion: bool,
                        with_luminosity_function: bool,
                        with_velocities_vs_magnitude: bool,
                        with_velocity_clouds: bool,
                        heatmaps_axes: str,
                        with_toomre_diagram: bool,
                        with_ugriz_diagrams: bool
                        ) -> List[InstrumentedAttribute]:
    entities = []

    if filtration_method != 'raw':
        entities += [Star.distance,
                     Star.declination,
                     Star.u_velocity,
                     Star.v_velocity,
                     Star.w_velocity]
    if filtration_method == 'restricted':
        entities += [Star.b_abs_magnitude,
                     Star.v_abs_magnitude,
                     Star.r_abs_magnitude,
                     Star.i_abs_magnitude,
                     Star.proper_motion]

    if nullify_radial_velocity:
        entities += [Star.galactic_longitude,
                     Star.galactic_latitude,
                     Star.proper_motion_component_l,
                     Star.proper_motion_component_b,
                     Star.distance]

    if lepine_criterion:
        entities += [Star.right_ascension,
                     Star.declination,
                     Star.distance]

    if with_luminosity_function:
        entities += [Star.luminosity]

    if with_velocities_vs_magnitude:
        entities += [Star.luminosity,
                     Star.u_velocity,
                     Star.v_velocity,
                     Star.w_velocity]

    if heatmaps_axes or with_velocity_clouds:
        entities += [Star.u_velocity,
                     Star.v_velocity,
                     Star.w_velocity]

    if with_toomre_diagram:
        entities += [Star.u_velocity,
                     Star.v_velocity,
                     Star.w_velocity,
                     Star.spectral_type]

    if with_ugriz_diagrams:
        entities += [Star.b_abs_magnitude,
                     Star.v_abs_magnitude,
                     Star.r_abs_magnitude,
                     Star.i_abs_magnitude,
                     Star.spectral_type]

    if entities:
        entities += [Star.id]

    return entities

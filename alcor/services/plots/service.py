import uuid
from collections import Counter
from functools import partial
from itertools import filterfalse

from sqlalchemy.orm.session import Session

from alcor.models import eliminations
from alcor.models.star import set_radial_velocity_to_zero
from alcor.services.data_access import fetch_group_stars
from alcor.services.stars_group import elimination
from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds,
               heatmaps,
               toomre_diagram,
               ugriz_diagrams)


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
         session: Session) -> None:
    if any((with_luminosity_function,
            with_velocity_clouds,
            with_velocities_vs_magnitude)):
        stars = fetch_group_stars(group_id=group_id,
                                  session=session)

        stars_count = len(stars)
        eliminations_counter = Counter()

        if filtration_method in {'restricted', 'full'}:
            is_eliminated = partial(elimination.check,
                                    eliminations_counter=eliminations_counter,
                                    filtration_method=filtration_method)
            stars = list(filterfalse(is_eliminated, stars))

        # TODO: figure out what to do when there is no group_id
        counter = eliminations.StarsCounter(group_id=group_id,
                                            raw=stars_count,
                                            **eliminations_counter)
        session.add(counter)
        session.commit()

        if nullify_radial_velocity:
            stars = list(map(set_radial_velocity_to_zero, stars))

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

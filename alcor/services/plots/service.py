import uuid
from typing import Optional

from sqlalchemy.orm.session import Session

from alcor.models import Star
from alcor.models.star import set_radial_velocity_to_zero
from alcor.services.data_access import (fetch_all,
                                        fetch_group_stars)
from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds,
               heatmaps,
               toomre_diagram,
               ugriz_diagrams)


def draw(group_id: Optional[uuid.UUID],
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
        if group_id:
            stars = fetch_group_stars(group_id=group_id,
                                      session=session)
        else:
            stars = fetch_all(Star,
                              session=session)
        if nullify_radial_velocity:
            stars = list(map(set_radial_velocity_to_zero, stars))

    if with_luminosity_function:
        luminosity_function.plot(stars)

    if with_velocities_vs_magnitude:
        if lepine_criterion:
            velocities_vs_magnitude.plot_lepine_case(stars)
        else:
            velocities_vs_magnitude.plot(stars)

    if with_velocity_clouds:
        if lepine_criterion:
            velocity_clouds.plot_lepine_case(stars)
        else:
            velocity_clouds.plot(stars)

    if heatmaps_axes:
        heatmaps.plot(group_id=group_id,
                      session=session,
                      axes=heatmaps_axes)

    if with_toomre_diagram:
        toomre_diagram.plot(session=session)

    if with_ugriz_diagrams:
        ugriz_diagrams.plot(session=session)

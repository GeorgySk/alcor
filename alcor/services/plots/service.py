import uuid
from typing import Optional

from sqlalchemy.orm.session import Session

from alcor.models.star import Star
from alcor.services.data_access import (fetch_all,
                                        fetch_group_stars)
from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds,
               heatmaps,
               toomre_diagram,
               ugriz_diagrams)


def draw(group_id: Optional[uuid.UUID],
         with_luminosity_function: bool,
         with_velocities_vs_magnitude: bool,
         with_velocity_clouds: bool,
         lepine_criterion: bool,
         heatmaps_axes: str,
         with_toomre_diagram: bool,
         with_ugriz_diagrams: bool,
         session: Session) -> None:
    if any((with_luminosity_function, with_velocity_clouds)):
        if group_id:
            stars = fetch_group_stars(group_id=group_id,
                                      session=session)
        else:
            stars = fetch_all(Star,
                              session=session)

    if with_luminosity_function:
        luminosity_function.plot(stars)

    if with_velocities_vs_magnitude:
        if lepine_criterion:
            velocities_vs_magnitude.plot_lepine_case(session=session)
        else:
            velocities_vs_magnitude.plot(session=session)

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

from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds,
               heatmaps)

from cassandra.cluster import Session


def draw_plots(with_luminosity_function: bool,
               with_velocities_vs_magnitude: bool,
               with_velocity_clouds: bool,
               lepine_criterion: bool,
               heatmaps_axes: str,
               session: Session) -> None:
    if with_luminosity_function:
        luminosity_function.plot(session=session)

    if with_velocities_vs_magnitude:
        if lepine_criterion:
            velocities_vs_magnitude.plot_lepine_case(session=session)
        else:
            velocities_vs_magnitude.plot(session=session)

    if with_velocity_clouds:
        if lepine_criterion:
            velocity_clouds.plot_lepine_case(session=session)
        else:
            velocity_clouds.plot(session=session)

    if heatmaps_axes:
        heatmaps.plot(session=session,
                      axes=heatmaps_axes)

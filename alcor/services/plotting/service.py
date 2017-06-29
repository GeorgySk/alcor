from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds)

from cassandra.cluster import Session


def draw_plots(luminosity_function_flag: bool,
               velocities_vs_magnitude_flag: bool,
               velocity_clouds_flag: bool,
               lepine_criterion: bool,
               session: Session) -> None:
    if luminosity_function_flag:
        luminosity_function.plot(session=session)

    if velocities_vs_magnitude_flag:
        if lepine_criterion:
            velocities_vs_magnitude.plot_lepine_case(session=session)
        else:
            velocities_vs_magnitude.plot(session=session)

    if velocity_clouds_flag:
        if lepine_criterion:
            velocity_clouds.plot_lepine_case(session=session)
        else:
            velocity_clouds.plot(session=session)

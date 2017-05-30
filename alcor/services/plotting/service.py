from . import (luminosity_function,
               velocities_vs_magnitude,
               velocity_clouds)


def draw_plots(luminosity_function_flag: bool,
               velocities_vs_magnitude_flag: bool,
               velocity_clouds_flag: bool,
               lepine_criterion: bool) -> None:
    if luminosity_function_flag:
        luminosity_function.plot()

    if velocities_vs_magnitude_flag and not lepine_criterion:
        velocities_vs_magnitude.plot()
    if velocities_vs_magnitude_flag and lepine_criterion:
        velocities_vs_magnitude.plot_lepine_case()

    if velocity_clouds_flag and not lepine_criterion:
        velocity_clouds.plot()
    if velocity_clouds_flag and lepine_criterion:
        velocity_clouds.plot_lepine_case()

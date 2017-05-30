from . import luminosity_function
from . import velocities_vs_magnitude


def make_plots(luminosity_function_flag: bool,
               velocities_vs_magnitude_flag: bool,
               lepine_criterion: bool) -> None:
    if luminosity_function_flag:
        luminosity_function.plot()

    if velocities_vs_magnitude_flag and not lepine_criterion:
        velocities_vs_magnitude.plot()
    if velocities_vs_magnitude_flag and lepine_criterion:
        velocities_vs_magnitude.plot_lepine_case()

from .luminosity_function import plot_luminosity_function
from .velocities_vs_magnitude import (
    plot_velocities_vs_magnitude,
    plot_velocities_vs_magnitude_lepine_case)


def make_plots(luminosity_function: bool,
               velocities_vs_magnitude: bool,
               lepine_criterion: bool) -> None:
    if luminosity_function:
        plot_luminosity_function()

    if velocities_vs_magnitude and not lepine_criterion:
        plot_velocities_vs_magnitude()
    if velocities_vs_magnitude and lepine_criterion:
        plot_velocities_vs_magnitude_lepine_case()

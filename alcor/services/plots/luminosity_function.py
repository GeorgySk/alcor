from math import (log10,
                  sqrt)
from typing import Tuple

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
from decimal import Decimal

matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from sqlalchemy.orm.session import Session

from alcor.models.luminosity_function import Point
from alcor.services.common import FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
from alcor.services.data_access import fetch_all

# Observational LF of 40pc sample from Althaus
OBSERVATIONAL_AVG_BIN_MAGNITUDES = np.arange(7.75, 17.25, 0.5)
OBSERVATIONAL_STARS_COUNTS = [3, 4, 5, 7, 12, 17, 17, 12, 20, 19, 37, 42, 52,
                              72, 96, 62, 20, 3, 1]
OBSERVATIONAL_STARS_COUNTS_LOGARITHMS = [
    log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
    for stars_count in OBSERVATIONAL_STARS_COUNTS]
OBSERVATIONAL_UPPER_ERRORBARS = [
    log10((stars_count + sqrt(stars_count))
          / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
    - log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
    for stars_count in OBSERVATIONAL_STARS_COUNTS]
EPSILON = 1e-6
OBSERVATIONAL_LOWER_ERRORBARS = [
    log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
    - log10((stars_count - sqrt(stars_count) + EPSILON)
            / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
    for stars_count in OBSERVATIONAL_STARS_COUNTS]
OBSERVATIONAL_ASYMMETRIC_ERRORBARS = [OBSERVATIONAL_LOWER_ERRORBARS,
                                      OBSERVATIONAL_UPPER_ERRORBARS]


def plot(session: Session,
         filename: str = 'luminosity_function.ps',
         figure_size: Tuple[float, float] = (7, 7),
         ratio: float = 10 / 13,
         xlabel: str = '$M_{bol}$',
         ylabel: str = '$\log N (pc^{-3}M_{bol}^{-1})$',
         xlimits: Tuple[float, float] = (7, 19),
         ylimits: Tuple[float, float] = (-6, -2),
         line_color: str = 'k',
         marker: str = 's',
         capsize: float = 5,
         observational_line_color: str = 'r') -> None:
    # TODO: Implement other fetching functions
    graph_points = fetch_all(Point,
                             session=session)

    avg_bin_magnitudes = [graph_point.avg_bin_magnitude
                          for graph_point in graph_points]
    stars_count_logarithms = [graph_point.stars_count_logarithm
                              for graph_point in graph_points]

    normalized_stars_counts = [
        FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
        * 10. ** float(stars_count_logarithm)
        for stars_count_logarithm in stars_count_logarithms]
    upper_errorbars = [
        log10((stars_count + sqrt(stars_count))
              / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
        - log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
        for stars_count in normalized_stars_counts]

    lower_errorbars = []
    for stars_count in normalized_stars_counts:
        try:
            lower_errorbars.append(
                    log10(stars_count
                          / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                    - log10((stars_count - sqrt(stars_count) + EPSILON)
                            / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))
        except ValueError:
            # Some number so that errorbar would go below the plot
            lower_errorbars.append(5.0)

    upper_errorbars = [Decimal(errorbar) for errorbar in upper_errorbars]
    lower_errorbars = [Decimal(errorbar) for errorbar in lower_errorbars]

    (avg_bin_magnitudes,
     stars_count_logarithms,
     upper_errorbars,
     lower_errorbars) = zip(*sorted(zip(avg_bin_magnitudes,
                                        stars_count_logarithms,
                                        upper_errorbars,
                                        lower_errorbars)))

    asymmetric_errorbars = [lower_errorbars,
                            upper_errorbars]

    figure, subplot = plt.subplots(figsize=figure_size)

    subplot.set(xlabel=xlabel,
                ylabel=ylabel,
                xlim=xlimits,
                ylim=ylimits)

    subplot.errorbar(x=avg_bin_magnitudes,
                     y=stars_count_logarithms,
                     yerr=asymmetric_errorbars,
                     marker=marker,
                     color=line_color,
                     capsize=capsize,
                     zorder=2)

    subplot.errorbar(x=OBSERVATIONAL_AVG_BIN_MAGNITUDES,
                     y=OBSERVATIONAL_STARS_COUNTS_LOGARITHMS,
                     yerr=OBSERVATIONAL_ASYMMETRIC_ERRORBARS,
                     marker=marker,
                     color=observational_line_color,
                     capsize=capsize,
                     zorder=1)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    plt.savefig(filename)

from math import (log10,
                  sqrt)

from sqlalchemy.orm.session import Session
import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
from alcor.services.data_access import fetch_all_graph_points

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from alcor.services.restrictions import FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME

FILENAME = 'luminosity_function.ps'

FIGURE_SIZE = (7, 7)
DESIRED_DIMENSIONS_RATIO = 10 / 13

X_LABEL = '$M_{bol}$'
Y_LABEL = '$\log N (pc^{-3}M_{bol}^{-1})$'

X_LIMITS = [7, 19]
Y_LIMITS = [-6, -2]

LINE_COLOR = 'k'
MARKER = 's'
CAP_SIZE = 5

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

OBSERVATIONAL_LINE_COLOR = 'r'


def plot(session: Session) -> None:
    # TODO: Implement other fetching functions
    graph_points = fetch_all_graph_points(session=session)

    avg_bin_magnitudes = [graph_point.avg_bin_magnitude
                          for graph_point in graph_points]
    stars_count_logarithms = [graph_point.stars_count_logarithm
                              for graph_point in graph_points]
    upper_errorbars = [graph_point.upper_error_bar
                       for graph_point in graph_points]
    lower_errorbars = [graph_point.lower_error_bar
                       for graph_point in graph_points]

    (avg_bin_magnitudes,
     stars_count_logarithms,
     upper_errorbars,
     lower_errorbars) = zip(*sorted(zip(avg_bin_magnitudes,
                                        stars_count_logarithms,
                                        upper_errorbars,
                                        lower_errorbars)))

    asymmetric_errorbars = [lower_errorbars,
                            upper_errorbars]

    figure, subplot = plt.subplots(figsize=FIGURE_SIZE)

    subplot.set(xlabel=X_LABEL,
                ylabel=Y_LABEL,
                xlim=X_LIMITS,
                ylim=Y_LIMITS)

    subplot.errorbar(x=avg_bin_magnitudes,
                     y=stars_count_logarithms,
                     yerr=asymmetric_errorbars,
                     marker=MARKER,
                     color=LINE_COLOR,
                     capsize=CAP_SIZE,
                     zorder=2)

    subplot.errorbar(x=OBSERVATIONAL_AVG_BIN_MAGNITUDES,
                     y=OBSERVATIONAL_STARS_COUNTS_LOGARITHMS,
                     yerr=OBSERVATIONAL_ASYMMETRIC_ERRORBARS,
                     marker=MARKER,
                     color=OBSERVATIONAL_LINE_COLOR,
                     capsize=CAP_SIZE,
                     zorder=1)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot.get_data_ratio())

    plt.savefig(FILENAME)

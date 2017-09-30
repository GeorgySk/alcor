from math import (log10,
                  sqrt)
from typing import (Iterable,
                    Iterator,
                    Tuple,
                    List)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
from alcor.types import StarsBinsType

matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from alcor.models import Star
from alcor.models.luminosity_function import Point
from alcor.services.common import (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME,
                                   stars_packer,
                                   star_bolometric_indexer,
                                   bolometric_indexer)

MIN_BOLOMETRIC_MAGNITUDE = 6.
STARS_BIN_SIZE = 0.5
MAX_BOLOMETRIC_MAGNITUDE = 21.
# More info at (The white dwarf population within 40 pc of the Sun,
# Torres et al., 2016)
SMALLEST_ERROR_BARS_BINS_INDEXES = {15, 16, 17}
OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT = 220
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
bolometric_index = bolometric_indexer(min_magnitude=MIN_BOLOMETRIC_MAGNITUDE,
                                      stars_bin_size=STARS_BIN_SIZE)
STARS_BINS_COUNT = bolometric_index(MAX_BOLOMETRIC_MAGNITUDE)
star_bolometric_index = star_bolometric_indexer(bolometric_index)
pack_stars = stars_packer(stars_bins_count=STARS_BINS_COUNT,
                          star_bolometric_index=star_bolometric_index)


def plot(stars: List[Star],
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
    stars_bins = pack_stars(stars)
    normalization_factor = stars_bins_normalization_factor(stars_bins)
    stars_counts = map(len, stars_bins)
    points = list(graph_points(stars_counts=stars_counts,
                               normalization_factor=normalization_factor))

    avg_bin_magnitudes = [graph_point.avg_bin_magnitude
                          for graph_point in points]
    stars_count_logarithms = [graph_point.stars_count_logarithm
                              for graph_point in points]

    # As simulated number of white dwarfs (WDs) can be much greater than the
    # number of WDs for which observational luminosity function (LF) is
    # plotted, we normalize the LF based on synthetic(simulated) WDs
    # to the observational one (move to the same location on the plot).
    # We do it in 'processing' module. TODO: move it here
    # Another thing is that error bars become very small for synthetic LF
    # and here we 1) get number of stars for which we would have the
    # corresponding normalized logN and 2) calculate new not squeezed error
    # bars corresponding to this number of stars
    normalized_stars_counts = [
        FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME * 10. ** stars_count_logarithm
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


def graph_points(*,
                 stars_counts: Iterable[int],
                 normalization_factor: float) -> Iterator[Point]:
    for index, stars_count in enumerate(stars_counts):
        if not stars_count:
            continue

        avg_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                             + STARS_BIN_SIZE * (index + 0.5))
        logarithm = to_logarithm(stars_count=stars_count,
                                 normalization_factor=normalization_factor)
        upper_error_bar = to_upper_error_bar(stars_count)
        lower_error_bar = to_lower_error_bar(stars_count)
        yield Point(avg_bin_magnitude=avg_bin_magnitude,
                    stars_count_logarithm=logarithm,
                    upper_error_bar=upper_error_bar,
                    lower_error_bar=lower_error_bar)


def stars_bins_normalization_factor(
        stars_bins: StarsBinsType,
        stars_bins_indexes: Iterable[int] = SMALLEST_ERROR_BARS_BINS_INDEXES
) -> float:
    smallest_error_bars_bins_lengths_sum = (
        sum(len(stars_bins[index])
            for index in stars_bins_indexes))
    return (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
            * smallest_error_bars_bins_lengths_sum
            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)


def to_logarithm(stars_count: int,
                 normalization_factor: float) -> float:
    try:
        return log10(stars_count / normalization_factor)
    except ValueError:
        return 0.


def to_upper_error_bar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    else:
        return (log10((stars_count + sqrt(stars_count))
                      / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))


def to_lower_error_bar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    elif stars_count == 1:
        return 5.  # Some number so that error bar would go below the plot
    else:
        return (log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10((stars_count - sqrt(stars_count))
                        / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))

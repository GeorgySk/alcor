from functools import partial
from typing import (Callable,
                    Tuple,
                    List)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

from alcor.models import Star

OBSERVATIONAL_STARS_COUNTS = np.array(
        [3, 4, 5, 7, 12, 17, 17, 12, 20, 19, 37, 42, 52, 72, 96, 62, 20, 3, 1])

# TODO: replace by # 2 * pi * radius ** 3 / 3
FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME = 134041.29
nan_array = partial(np.full,
                    fill_value=np.nan)


def bolometric_indexer(*,
                       min_magnitude: float,
                       stars_bin_size: float) -> Callable[[np.ndarray],
                                                          np.ndarray]:
    def bolometric_index(magnitude: np.ndarray) -> np.ndarray:
        magnitude_amplitude = magnitude - min_magnitude
        return np.floor(magnitude_amplitude / stars_bin_size).astype(np.int32)

    return bolometric_index


# TODO: there is no need to have stars with all the fields
def plot(stars: List[Star],
         min_bolometric_magnitude: float = 6.,
         max_bolometric_magnitude: float = 21.,
         bin_size: float = 0.5,
         min_observed_magnitude: float = 7.75,
         observed_stars_counts: np.ndarray = OBSERVATIONAL_STARS_COUNTS,
         # We choose these bins because they have many objects and don't lie in
         # problematic regions
         trusted_bins: frozenset = {15, 16, 17},
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
    bolometric_index = bolometric_indexer(
        min_magnitude=min_bolometric_magnitude,
        stars_bin_size=bin_size)

    stars_bins_count = np.asscalar(bolometric_index(max_bolometric_magnitude))

    # Aligning observed stars counts with the scale defined by min and max
    # bolometric magnitudes
    initial_index = np.asscalar(bolometric_index(min_observed_magnitude))
    observed_stars_counts = np.insert(arr=observed_stars_counts,
                                      obj=0,
                                      values=np.zeros(shape=initial_index,
                                                      dtype=np.int32))
    observed_stars_counts = np.append(
            arr=observed_stars_counts,
            values=np.zeros(
                    shape=stars_bins_count - observed_stars_counts.size,
                    dtype=np.int32))

    observational_luminosity_function = luminosity_function(
            max_bolometric_magnitude=max_bolometric_magnitude,
            min_bolometric_magnitude=min_bolometric_magnitude,
            bin_size=bin_size,
            stars_bins_count=stars_bins_count,
            stars_counts=observed_stars_counts)

    magnitudes = np.array([star.bolometric_magnitude
                           for star in stars])
    bins_indexes = pd.Series(bolometric_index(magnitudes))

    bins_counters = np.zeros(shape=stars_bins_count,
                             dtype=np.int32)
    # More info on np.ndarray iterations at:
    # docs.scipy.org/doc/numpy-1.13.0/reference/arrays.nditer.html#modifying-array-values
    # More info on indexation at:
    # docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.ndindex.html
    for index, bin_counter in zip(np.ndindex(bins_counters.shape),
                                  np.nditer(op=bins_counters,
                                            op_flags=['readwrite'])):
        bin_counter[...] = bins_indexes[bins_indexes == index].count()

    # As simulated number of white dwarfs (WDs) can be much greater than the
    # number of WDs for which observational luminosity function (LF) is
    # plotted, we normalize the LF based on synthetic(simulated) WDs
    # to the observational one (move to the same location on the plot).
    actual_stars_counts = pd.Series(bins_counters)

    observed_stars_counts = pd.Series(observed_stars_counts)
    trusted_bins_observed_stars_count = (
        observed_stars_counts[
            observed_stars_counts.index.isin(trusted_bins)].sum())
    trusted_bins_simulated_stars_count = (
        actual_stars_counts[
            actual_stars_counts.index.isin(trusted_bins)].sum())

    normalized_stars_counts = (actual_stars_counts
                               * trusted_bins_observed_stars_count
                               / trusted_bins_simulated_stars_count)

    synthetic_luminosity_function = luminosity_function(
            max_bolometric_magnitude=max_bolometric_magnitude,
            min_bolometric_magnitude=min_bolometric_magnitude,
            bin_size=bin_size,
            stars_bins_count=stars_bins_count,
            stars_counts=normalized_stars_counts)

    figure, subplot = plt.subplots(figsize=figure_size)

    subplot.set(xlabel=xlabel,
                ylabel=ylabel,
                xlim=xlimits,
                ylim=ylimits)

    subplot.errorbar(x=synthetic_luminosity_function['magnitude'],
                     y=synthetic_luminosity_function['log_stars_count'],
                     yerr=[synthetic_luminosity_function['lower_errorbar'],
                           synthetic_luminosity_function['upper_errorbar']],
                     marker=marker,
                     color=line_color,
                     capsize=capsize,
                     zorder=2)

    subplot.errorbar(x=observational_luminosity_function['magnitude'],
                     y=observational_luminosity_function['log_stars_count'],
                     yerr=[
                         observational_luminosity_function['lower_errorbar'],
                         observational_luminosity_function['upper_errorbar']],
                     marker=marker,
                     color=observational_line_color,
                     capsize=capsize,
                     zorder=1)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    plt.savefig(filename)


def luminosity_function(
        min_bolometric_magnitude: float,
        max_bolometric_magnitude: float,
        bin_size: float,
        stars_bins_count: int,
        stars_counts: np.ndarray,
        max_errorbar_len: float = 6.
) -> DataFrame:
    luminosity_function_template = dict(
            magnitude=np.arange(min_bolometric_magnitude + bin_size / 2,
                                max_bolometric_magnitude,
                                bin_size),
            log_stars_count=nan_array(stars_bins_count),
            upper_errorbar=nan_array(stars_bins_count),
            lower_errorbar=nan_array(stars_bins_count))

    res = pd.DataFrame(data=luminosity_function_template)

    res['log_stars_count'] = np.log10(
            stars_counts / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
    res['upper_errorbar'] = np.log10(1. + 1. / np.sqrt(stars_counts))
    res['lower_errorbar'] = -np.log10(1. - 1. / np.sqrt(stars_counts))

    res.replace(to_replace=[np.inf, -np.inf],
                value=np.nan,
                inplace=True)

    # Replacing erroneous NaN values for errorbars
    missing_values_rows = res.loc[
        (res.isnull().any(axis=1))
        & (res['log_stars_count'].notnull()), :].copy()
    missing_values_rows.fillna(value=max_errorbar_len,
                               inplace=True)
    res.loc[(res.isnull().any(axis=1))
            & (res['log_stars_count'].notnull()), :] = missing_values_rows

    return res

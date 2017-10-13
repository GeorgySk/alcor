from functools import partial
from typing import Tuple

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from .utils import (bolometric_indexer,
                    bolometric_magnitude,
                    nan_array)

OBSERVATIONAL_STARS_COUNTS = np.array(
        [3, 4, 5, 7, 12, 17, 17, 12, 20, 19, 37, 42, 52, 72, 96, 62, 20, 3, 1])

# TODO: replace by # 2 * pi * radius ** 3 / 3
FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME = 134041.29


def plot(stars: pd.DataFrame,
         *,
         min_bolometric_magnitude: float = 6.,
         max_bolometric_magnitude: float = 21.,
         bin_size: float = 0.5,
         min_observed_magnitude: float = 7.75,
         observed_stars_counts: np.ndarray = OBSERVATIONAL_STARS_COUNTS,
         # We choose these bins because they have many objects
         # and don't lie in problematic regions
         trusted_bins: frozenset = frozenset([15, 16, 17]),
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

    # Aligning observed stars counts with the scale
    # defined by min and max bolometric magnitudes
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

    magnitudes = bolometric_magnitude(luminosities=stars['luminosity'])
    bins_indexes = pd.Series(bolometric_index(magnitudes))

    actual_stars_counts = bins_indexes.value_counts()

    observed_stars_counts = pd.Series(observed_stars_counts)

    normalized_stars_counts = actual_stars_counts * normalization_factor(
            actual_stars_counts=actual_stars_counts,
            observed_stars_counts=observed_stars_counts,
            trusted_bins=trusted_bins)

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

    draw_errorbar = partial(subplot.errorbar,
                            marker=marker,
                            capsize=capsize)

    draw_errorbar(x=synthetic_luminosity_function['magnitude'],
                  y=synthetic_luminosity_function['log_stars_count'],
                  yerr=[synthetic_luminosity_function['lower_errorbar'],
                        synthetic_luminosity_function['upper_errorbar']],
                  color=line_color,
                  zorder=2)

    draw_errorbar(x=observational_luminosity_function['magnitude'],
                  y=observational_luminosity_function['log_stars_count'],
                  yerr=[observational_luminosity_function['lower_errorbar'],
                        observational_luminosity_function['upper_errorbar']],
                  color=observational_line_color,
                  zorder=1)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    plt.savefig(filename)


def luminosity_function(*,
                        min_bolometric_magnitude: float,
                        max_bolometric_magnitude: float,
                        bin_size: float,
                        stars_bins_count: int,
                        stars_counts: np.ndarray,
                        max_errorbar_len: float = 6.) -> pd.DataFrame:
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

    stars_counts_sqrt = np.sqrt(stars_counts)
    inverse_stars_counts_sqrt = 1. / stars_counts_sqrt
    res['upper_errorbar'] = np.log10(1. + inverse_stars_counts_sqrt)
    res['lower_errorbar'] = -np.log10(1. - inverse_stars_counts_sqrt)

    res.replace(to_replace=[np.inf, -np.inf],
                value=np.nan,
                inplace=True)

    return replace_nans(df=res,
                        replacement=max_errorbar_len)


def trusted_bins_stars_count(stars_counts: pd.Series,
                             *,
                             trusted_bins: frozenset) -> float:
    return stars_counts[stars_counts.index.isin(trusted_bins)].sum()


def normalization_factor(*,
                         trusted_bins: frozenset,
                         observed_stars_counts: pd.Series,
                         actual_stars_counts: pd.Series) -> float:
    trusted_bins_stars_counter = partial(trusted_bins_stars_count,
                                         trusted_bins=trusted_bins)

    return (trusted_bins_stars_counter(observed_stars_counts)
            / trusted_bins_stars_counter(actual_stars_counts))


def replace_nans(df: pd.DataFrame,
                 replacement: float) -> pd.DataFrame:
    result = pd.DataFrame(data=df,
                          copy=True)

    missing_values_rows_mask = result.isnull().any(axis=1)
    notnull_log_stars_count_rows_mask = result['log_stars_count'].notnull()
    missing_values_rows_mask &= notnull_log_stars_count_rows_mask

    missing_values_rows = result.loc[missing_values_rows_mask, :].copy()

    missing_values_rows.fillna(value=replacement,
                               inplace=True)
    result.loc[missing_values_rows_mask, :] = missing_values_rows

    return result

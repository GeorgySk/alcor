from typing import (Callable,
                    Dict,
                    Tuple)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

from alcor.utils import zip_mappings
from .utils import (bolometric_indexer,
                    bolometric_magnitude,
                    nan_array,
                    to_cartesian_from_equatorial)


def plot(stars: pd.DataFrame,
         *,
         min_bolometric_magnitude: float = 6.,
         max_bolometric_magnitude: float = 30.,
         bin_size: float = 0.5,
         figure_size: Tuple[float, float] = (10, 12),
         filename: str = 'velocities_vs_magnitude.ps',
         u_label: str = '$U_{LSR}(km/s)$',
         v_label: str = '$V_{LSR}(km/s)$',
         w_label: str = '$W_{LSR}(km/s)$',
         magnitude_label: str = '$M_{bol}$') -> None:
    labels = dict(u_velocity=u_label,
                  v_velocity=v_label,
                  w_velocity=w_label)

    bolometric_index = bolometric_indexer(
            min_magnitude=min_bolometric_magnitude,
            stars_bin_size=bin_size)

    bins_by_velocities = empty_bins_by_velocities(
            min_bolometric_magnitude=min_bolometric_magnitude,
            max_bolometric_magnitude=max_bolometric_magnitude,
            bin_size=bin_size,
            bolometric_index=bolometric_index)

    magnitudes = bolometric_magnitude(luminosities=stars['luminosity'])
    bins_indexes = pd.Series(bolometric_index(magnitudes))

    figure, subplots = plt.subplots(nrows=3,
                                    figsize=figure_size)
    subplots = dict(u_velocity=subplots[0],
                    v_velocity=subplots[1],
                    w_velocity=subplots[2])

    for velocity, (bins,
                   subplot,
                   label) in zip_mappings(bins_by_velocities,
                                          subplots,
                                          labels):
        bins = fill_bins(bins,
                         stars=stars,
                         bins_indexes=bins_indexes,
                         velocity=velocity)

        draw_subplot(subplot=subplot,
                     ylabel=label,
                     x_line=bins['magnitude'],
                     y_line=bins['avg_velocity'],
                     yerr=bins['velocity_std'],
                     x_scatter=magnitudes,
                     y_scatter=stars[velocity])

    # Removing unnecessary x-labels for top and middle subplots
    subplots['u_velocity'].set_xticklabels([])
    subplots['v_velocity'].set_xticklabels([])

    subplots['w_velocity'].set_xlabel(magnitude_label)

    # TODO: delete overlapping y-labels
    figure.subplots_adjust(hspace=0)

    plt.savefig(filename)


def plot_lepine_case(stars: pd.DataFrame,
                     *,
                     min_bolometric_magnitude: float = 6.,
                     max_bolometric_magnitude: float = 30.,
                     bin_size: float = 0.5,
                     figure_size: Tuple[float, float] = (10, 12),
                     filename: str = 'velocities_vs_magnitude.ps',
                     u_label: str = '$U_{LSR}(km/s)$',
                     v_label: str = '$V_{LSR}(km/s)$',
                     w_label: str = '$W_{LSR}(km/s)$',
                     magnitude_label: str = '$M_{bol}$') -> None:
    labels = dict(u_velocity=u_label,
                  v_velocity=v_label,
                  w_velocity=w_label)

    bolometric_index = bolometric_indexer(
            min_magnitude=min_bolometric_magnitude,
            stars_bin_size=bin_size)

    bins_by_velocities = empty_bins_by_velocities(
            min_bolometric_magnitude=min_bolometric_magnitude,
            max_bolometric_magnitude=max_bolometric_magnitude,
            bin_size=bin_size,
            bolometric_index=bolometric_index)

    stars_by_velocities = split_stars_by_velocities(stars)

    figure, subplots = plt.subplots(nrows=3,
                                    figsize=figure_size)
    subplots = dict(u_velocity=subplots[0],
                    v_velocity=subplots[1],
                    w_velocity=subplots[2])

    for velocity, (bins,
                   stars,
                   subplot,
                   label) in zip_mappings(bins_by_velocities,
                                          stars_by_velocities,
                                          subplots,
                                          labels):
        magnitudes = bolometric_magnitude(luminosities=stars['luminosity'])
        bins_indexes = pd.Series(bolometric_index(magnitudes))

        bins = fill_bins(bins,
                         stars=stars,
                         bins_indexes=bins_indexes,
                         velocity=velocity)

        draw_subplot(subplot=subplot,
                     ylabel=label,
                     x_line=bins['magnitude'],
                     y_line=bins['avg_velocity'],
                     yerr=bins['velocity_std'],
                     x_scatter=magnitudes,
                     y_scatter=stars[velocity])

    # Removing unnecessary x-labels for top and middle subplots
    subplots['u_velocity'].set_xticklabels([])
    subplots['v_velocity'].set_xticklabels([])

    subplots['w_velocity'].set_xlabel(magnitude_label)

    # TODO: delete overlapping y-labels
    figure.subplots_adjust(hspace=0)

    plt.savefig(filename)


def draw_subplot(*,
                 subplot: Axes,
                 xlabel: str = None,
                 ylabel: str,
                 xlim: Tuple[float, float] = (6, 19),
                 ylim: Tuple[float, float] = (-150, 150),
                 x_line: pd.Series,
                 y_line: pd.Series,
                 yerr: pd.Series,
                 marker: str = 's',
                 markersize: float = 3.,
                 line_color: str = 'k',
                 capsize: float = 5.,
                 linewidth: float = 1.,
                 x_scatter: pd.Series,
                 y_scatter: pd.Series,
                 scatter_color: str = 'gray',
                 scatter_point_size: float = 1.,
                 ratio: float = 7 / 13) -> None:
    subplot.set(xlabel=xlabel,
                ylabel=ylabel,
                xlim=xlim,
                ylim=ylim)
    subplot.errorbar(x=x_line,
                     y=y_line,
                     yerr=yerr,
                     marker=marker,
                     markersize=markersize,
                     color=line_color,
                     capsize=capsize,
                     linewidth=linewidth)
    subplot.scatter(x=x_scatter,
                    y=y_scatter,
                    color=scatter_color,
                    s=scatter_point_size)

    subplot.minorticks_on()
    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')
    subplot.set_aspect(ratio / subplot.get_data_ratio())


def empty_bins_by_velocities(min_bolometric_magnitude: float,
                             max_bolometric_magnitude: float,
                             bin_size: float,
                             bolometric_index: Callable
                             ) -> Dict[str, pd.DataFrame]:
    stars_bins_count = np.asscalar(bolometric_index(max_bolometric_magnitude))

    bins_template = dict(
            magnitude=np.arange(min_bolometric_magnitude + bin_size / 2,
                                max_bolometric_magnitude,
                                bin_size),
            avg_velocity=nan_array(stars_bins_count),
            velocity_std=nan_array(stars_bins_count))

    u_vs_magnitude_bins = pd.DataFrame(data=bins_template)
    v_vs_magnitude_bins = pd.DataFrame(data=bins_template)
    w_vs_magnitude_bins = pd.DataFrame(data=bins_template)

    return dict(u_velocity=u_vs_magnitude_bins,
                v_velocity=v_vs_magnitude_bins,
                w_velocity=w_vs_magnitude_bins)


def fill_bins(bins: pd.DataFrame,
              *,
              stars: pd.DataFrame,
              bins_indexes: pd.Series,
              velocity: str) -> pd.DataFrame:
    # TODO: check cases when there are no stars in a bin
    for index in range(bins.shape[0]):
        bin_corresponding_stars_mask = bins_indexes == index
        bin_corresponding_stars = stars[bin_corresponding_stars_mask]
        velocities = bin_corresponding_stars[velocity]

        bins.loc[index, 'avg_velocity'] = velocities.mean()
        # FIXME: one star in a bin returns nothing for std
        bins.loc[index, 'velocity_std'] = velocities.std()

    return bins


def split_stars_by_velocities(stars: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    x_coordinates, y_coordinates, z_coordinates = to_cartesian_from_equatorial(
            stars)

    highest_coordinates = np.maximum.reduce([np.abs(x_coordinates),
                                             np.abs(y_coordinates),
                                             np.abs(z_coordinates)])

    x_highest_coordinate_mask = highest_coordinates == x_coordinates
    y_highest_coordinate_mask = highest_coordinates == y_coordinates
    z_highest_coordinate_mask = highest_coordinates == z_coordinates

    u_vs_magnitude_stars = stars[y_highest_coordinate_mask
                                 | z_highest_coordinate_mask]
    v_vs_magnitude_stars = stars[x_highest_coordinate_mask
                                 | z_highest_coordinate_mask]
    w_vs_magnitude_stars = stars[x_highest_coordinate_mask
                                 | y_highest_coordinate_mask]
    return dict(u_velocity=u_vs_magnitude_stars,
                v_velocity=v_vs_magnitude_stars,
                w_velocity=w_vs_magnitude_stars)

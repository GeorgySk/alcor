from math import radians
from typing import Tuple

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

from .utils import (bolometric_indexer,
                    bolometric_magnitude,
                    nan_array)


DEC_GPOLE = radians(27.128336)
RA_GPOLE = radians(192.859508)
AUX_ANGLE = radians(122.932)


def plot(stars: pd.DataFrame,
         min_bolometric_magnitude: float = 6.,
         max_bolometric_magnitude: float = 30.,
         bin_size: float = 0.5,
         figure_size: Tuple[float, float] = (10, 12),
         filename: str = 'velocities_vs_magnitude.ps',
         u_label: str = '$U_{LSR}(km/s)$',
         v_label: str = '$V_{LSR}(km/s)$',
         w_label: str = '$W_{LSR}(km/s)$',
         magnitude_label: str = '$M_{bol}$') -> None:
    bolometric_index = bolometric_indexer(
            min_magnitude=min_bolometric_magnitude,
            stars_bin_size=bin_size)

    stars_bins_count = np.asscalar(bolometric_index(max_bolometric_magnitude))

    bins_template = dict(
            magnitude=np.arange(min_bolometric_magnitude + bin_size / 2,
                                max_bolometric_magnitude,
                                bin_size),
            avg_u_velocity=nan_array(stars_bins_count),
            avg_v_velocity=nan_array(stars_bins_count),
            avg_w_velocity=nan_array(stars_bins_count),
            u_velocity_std=nan_array(stars_bins_count),
            v_velocity_std=nan_array(stars_bins_count),
            w_velocity_std=nan_array(stars_bins_count))

    bins = pd.DataFrame(data=bins_template)

    magnitudes = bolometric_magnitude(luminosities=stars['luminosity'])
    bins_indexes = pd.Series(bolometric_index(magnitudes))

    # TODO: check cases when there are no stars in a bin
    for index in range(stars_bins_count):
        bins.loc[index, 'avg_u_velocity'] = stars[
            bins_indexes == index]['u_velocity'].mean()
        bins.loc[index, 'avg_v_velocity'] = stars[
            bins_indexes == index]['v_velocity'].mean()
        bins.loc[index, 'avg_w_velocity'] = stars[
            bins_indexes == index]['w_velocity'].mean()
        bins.loc[index, 'u_velocity_std'] = stars[
            bins_indexes == index]['u_velocity'].std()
        bins.loc[index, 'v_velocity_std'] = stars[
            bins_indexes == index]['v_velocity'].std()
        bins.loc[index, 'w_velocity_std'] = stars[
            bins_indexes == index]['w_velocity'].std()

    figure, (subplot_u,
             subplot_v,
             subplot_w) = plt.subplots(nrows=3,
                                       figsize=figure_size)

    draw_subplot(subplot=subplot_u,
                 ylabel=u_label,
                 x_line=bins['magnitude'],
                 y_line=bins['avg_u_velocity'],
                 yerr=bins['u_velocity_std'],
                 x_scatter=magnitudes,
                 y_scatter=stars['u_velocity'])
    draw_subplot(subplot=subplot_v,
                 ylabel=v_label,
                 x_line=bins['magnitude'],
                 y_line=bins['avg_v_velocity'],
                 yerr=bins['v_velocity_std'],
                 x_scatter=magnitudes,
                 y_scatter=stars['v_velocity'])
    draw_subplot(subplot=subplot_w,
                 xlabel=magnitude_label,
                 ylabel=w_label,
                 x_line=bins['magnitude'],
                 y_line=bins['avg_w_velocity'],
                 yerr=bins['w_velocity_std'],
                 x_scatter=magnitudes,
                 y_scatter=stars['w_velocity'])

    # Removing unnecessary x-labels for top and middle subplots
    subplot_u.set_xticklabels([])
    subplot_v.set_xticklabels([])

    # TODO: delete overlapping y-labels
    figure.subplots_adjust(hspace=0)

    plt.savefig(filename)


def plot_lepine_case(stars: pd.DataFrame,
                     min_bolometric_magnitude: float = 6.,
                     max_bolometric_magnitude: float = 30.,
                     bin_size: float = 0.5,
                     figure_size: Tuple[float, float] = (10, 12),
                     filename: str = 'velocities_vs_magnitude.ps',
                     u_label: str = '$U_{LSR}(km/s)$',
                     v_label: str = '$V_{LSR}(km/s)$',
                     w_label: str = '$W_{LSR}(km/s)$',
                     magnitude_label: str = '$M_{bol}$') -> None:
    bolometric_index = bolometric_indexer(
            min_magnitude=min_bolometric_magnitude,
            stars_bin_size=bin_size)

    stars_bins_count = np.asscalar(bolometric_index(max_bolometric_magnitude))

    bins_template = dict(
            magnitude=np.arange(min_bolometric_magnitude + bin_size / 2,
                                max_bolometric_magnitude,
                                bin_size),
            avg_velocity=nan_array(stars_bins_count),
            velocity_std=nan_array(stars_bins_count))

    u_vs_mag_bins = pd.DataFrame(data=bins_template)
    v_vs_mag_bins = pd.DataFrame(data=bins_template)
    w_vs_mag_bins = pd.DataFrame(data=bins_template)

    x_coordinates, y_coordinates, z_coordinates = to_cartesian_from_equatorial(
            stars)

    highest_coordinates = np.maximum.reduce([np.abs(x_coordinates),
                                             np.abs(y_coordinates),
                                             np.abs(z_coordinates)])

    u_vs_mag_stars = stars[(highest_coordinates == y_coordinates)
                           | (highest_coordinates == z_coordinates)]
    v_vs_mag_stars = stars[(highest_coordinates == x_coordinates)
                           | (highest_coordinates == z_coordinates)]
    w_vs_mag_stars = stars[(highest_coordinates == x_coordinates)
                           | (highest_coordinates == y_coordinates)]

    u_magnitudes = bolometric_magnitude(
            luminosity=u_vs_mag_stars['luminosity'])
    v_magnitudes = bolometric_magnitude(
            luminosity=v_vs_mag_stars['luminosity'])
    w_magnitudes = bolometric_magnitude(
            luminosity=w_vs_mag_stars['luminosity'])
    u_bins_indexes = pd.Series(bolometric_index(u_magnitudes))
    v_bins_indexes = pd.Series(bolometric_index(v_magnitudes))
    w_bins_indexes = pd.Series(bolometric_index(w_magnitudes))

    # TODO: check cases when there are no stars in a bin
    for index in range(stars_bins_count):
        u_vs_mag_bins.loc[index, 'avg_velocity'] = u_vs_mag_stars[
            u_bins_indexes == index]['u_velocity'].mean()
        v_vs_mag_bins.loc[index, 'avg_velocity'] = v_vs_mag_stars[
            v_bins_indexes == index]['v_velocity'].mean()
        w_vs_mag_bins.loc[index, 'avg_velocity'] = w_vs_mag_stars[
            w_bins_indexes == index]['w_velocity'].mean()
        u_vs_mag_bins.loc[index, 'velocity_std'] = u_vs_mag_stars[
            u_bins_indexes == index]['u_velocity'].std()
        v_vs_mag_bins.loc[index, 'velocity_std'] = v_vs_mag_stars[
            v_bins_indexes == index]['v_velocity'].std()
        w_vs_mag_bins.loc[index, 'velocity_std'] = w_vs_mag_stars[
            w_bins_indexes == index]['w_velocity'].std()

    figure, (subplot_u,
             subplot_v,
             subplot_w) = plt.subplots(nrows=3,
                                       figsize=figure_size)

    draw_subplot(subplot=subplot_u,
                 ylabel=u_label,
                 x_line=u_vs_mag_bins['magnitude'],
                 y_line=u_vs_mag_bins['avg_velocity'],
                 yerr=u_vs_mag_bins['velocity_std'],
                 x_scatter=u_magnitudes,
                 y_scatter=u_vs_mag_stars['u_velocity'])
    draw_subplot(subplot=subplot_v,
                 ylabel=v_label,
                 x_line=v_vs_mag_bins['magnitude'],
                 y_line=v_vs_mag_bins['avg_velocity'],
                 yerr=v_vs_mag_bins['velocity_std'],
                 x_scatter=v_magnitudes,
                 y_scatter=v_vs_mag_stars['v_velocity'])
    draw_subplot(subplot=subplot_w,
                 xlabel=magnitude_label,
                 ylabel=w_label,
                 x_line=w_vs_mag_bins['magnitude'],
                 y_line=w_vs_mag_bins['avg_velocity'],
                 yerr=w_vs_mag_bins['velocity_std'],
                 x_scatter=w_magnitudes,
                 y_scatter=w_vs_mag_stars['w_velocity'])

    # Removing unnecessary x-labels for top and middle subplots
    subplot_u.set_xticklabels([])
    subplot_v.set_xticklabels([])

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


def to_cartesian_from_equatorial(stars: pd.DataFrame) -> Tuple[pd.Series,
                                                               pd.Series,
                                                               pd.Series]:
    latitudes = (np.arcsin(np.cos(stars['declination']) * np.cos(DEC_GPOLE)
                           * np.cos(stars['right_ascension'] - RA_GPOLE)
                           + np.sin(stars['declination']) * np.sin(DEC_GPOLE)))
    x = np.sin(stars['declination']) - np.sin(latitudes) * np.sin(DEC_GPOLE)
    y = (np.cos(stars['declination'])
         * np.sin(stars['right_ascension'] - RA_GPOLE) * np.cos(DEC_GPOLE))
    longitudes = np.arctan(x / y) + AUX_ANGLE - np.pi / 2.
    longitudes[((x > 0.) & (y < 0.)) | ((x <= 0.) & (y <= 0.))] += np.pi

    x_coordinates = stars['distance'] * np.cos(latitudes) * np.cos(longitudes)
    y_coordinates = stars['distance'] * np.cos(latitudes) * np.sin(longitudes)
    z_coordinates = stars['distance'] * np.sin(latitudes)
    return x_coordinates, y_coordinates, z_coordinates

from typing import (Tuple,
                    List)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.axes import Axes

from alcor.models import Star
from alcor.services.velocities_vs_magnitudes import (lepine,
                                                     raw)


def plot(stars: List[Star],
         figure_size: Tuple[float, float] = (10, 12),
         filename: str = 'velocities_vs_magnitude.ps',
         u_label: str = '$U_{LSR}(km/s)$',
         v_label: str = '$V_{LSR}(km/s)$',
         w_label: str = '$W_{LSR}(km/s)$',
         magnitude_label: str = '$M_{bol}$') -> None:
    figure, (subplot_u,
             subplot_v,
             subplot_w) = plt.subplots(nrows=3,
                                       figsize=figure_size)

    bins = list(raw.bins(stars=stars))

    avg_bin_magnitudes = [float(stars_bin.avg_magnitude)
                          for stars_bin in bins]
    avg_u_velocities = [float(stars_bin.avg_u_velocity)
                        for stars_bin in bins]
    avg_velocities_v = [float(stars_bin.avg_v_velocity)
                        for stars_bin in bins]
    avg_velocities_w = [float(stars_bin.avg_w_velocity)
                        for stars_bin in bins]
    u_velocities_std = [float(stars_bin.u_velocity_std)
                        for stars_bin in bins]
    velocities_v_std = [float(stars_bin.v_velocity_std)
                        for stars_bin in bins]
    velocities_w_std = [float(stars_bin.w_velocity_std)
                        for stars_bin in bins]

    magnitudes = [float(star.bolometric_magnitude)
                  for star in stars]
    u_velocities = [float(star.u_velocity)
                    for star in stars]
    v_velocities = [float(star.v_velocity)
                    for star in stars]
    w_velocities = [float(star.w_velocity)
                    for star in stars]

    draw_subplot(subplot=subplot_u,
                 ylabel=u_label,
                 x_line=avg_bin_magnitudes,
                 y_line=avg_u_velocities,
                 yerr=u_velocities_std,
                 x_scatter=magnitudes,
                 y_scatter=u_velocities)
    draw_subplot(subplot=subplot_v,
                 ylabel=v_label,
                 x_line=avg_bin_magnitudes,
                 y_line=avg_velocities_v,
                 yerr=velocities_v_std,
                 x_scatter=magnitudes,
                 y_scatter=v_velocities)
    draw_subplot(subplot=subplot_w,
                 xlabel=magnitude_label,
                 ylabel=w_label,
                 x_line=avg_bin_magnitudes,
                 y_line=avg_velocities_w,
                 yerr=velocities_w_std,
                 x_scatter=magnitudes,
                 y_scatter=w_velocities)

    # Removing unnecessary x-labels for top and middle subplots
    subplot_u.set_xticklabels([])
    subplot_v.set_xticklabels([])

    # TODO: delete overlapping y-labels
    figure.subplots_adjust(hspace=0)

    plt.savefig(filename)


def plot_lepine_case(stars: List[Star],
                     figure_size: Tuple[float, float] = (10, 12),
                     filename: str = 'velocities_vs_magnitude.ps',
                     u_label: str = '$U_{LSR}(km/s)$',
                     v_label: str = '$V_{LSR}(km/s)$',
                     w_label: str = '$W_{LSR}(km/s)$',
                     magnitude_label: str = '$M_{bol}$') -> None:
    figure, (subplot_u,
             subplot_v,
             subplot_w) = plt.subplots(nrows=3,
                                       figsize=figure_size)

    u_vs_mag_bins, v_vs_mag_bins, w_vs_mag_bins = lepine.bins(stars=stars)

    u_bins_avg_magnitudes = [float(stars_bin.avg_magnitude)
                             for stars_bin in u_vs_mag_bins]
    avg_u_velocities = [float(stars_bin.avg_u_velocity)
                        for stars_bin in u_vs_mag_bins]
    u_velocities_std = [float(stars_bin.u_velocity_std)
                        for stars_bin in u_vs_mag_bins]
    v_bins_avg_magnitudes = [float(stars_bin.avg_magnitude)
                             for stars_bin in v_vs_mag_bins]
    avg_velocities_v = [float(stars_bin.avg_v_velocity)
                        for stars_bin in v_vs_mag_bins]
    velocities_v_std = [float(stars_bin.v_velocity_std)
                        for stars_bin in v_vs_mag_bins]
    w_bins_avg_magnitudes = [float(stars_bin.avg_magnitude)
                             for stars_bin in w_vs_mag_bins]
    avg_velocities_w = [float(stars_bin.avg_w_velocity)
                        for stars_bin in w_vs_mag_bins]
    velocities_w_std = [float(stars_bin.w_velocity_std)
                        for stars_bin in w_vs_mag_bins]

    u_vs_mag_cloud, v_vs_mag_cloud, w_vs_mag_cloud = lepine.clouds(stars=stars)

    u_magnitudes = [float(star.bolometric_magnitude)
                    for star in u_vs_mag_cloud]
    u_velocities = [float(star.u_velocity)
                    for star in u_vs_mag_cloud]
    v_magnitudes = [float(star.bolometric_magnitude)
                    for star in v_vs_mag_cloud]
    v_velocities = [float(star.v_velocity)
                    for star in v_vs_mag_cloud]
    w_magnitudes = [float(star.bolometric_magnitude)
                    for star in w_vs_mag_cloud]
    w_velocities = [float(star.w_velocity)
                    for star in w_vs_mag_cloud]

    draw_subplot(subplot=subplot_u,
                 ylabel=u_label,
                 x_line=u_bins_avg_magnitudes,
                 y_line=avg_u_velocities,
                 yerr=u_velocities_std,
                 x_scatter=u_magnitudes,
                 y_scatter=u_velocities)
    draw_subplot(subplot=subplot_v,
                 ylabel=v_label,
                 x_line=v_bins_avg_magnitudes,
                 y_line=avg_velocities_v,
                 yerr=velocities_v_std,
                 x_scatter=v_magnitudes,
                 y_scatter=v_velocities)
    draw_subplot(subplot=subplot_w,
                 xlabel=magnitude_label,
                 ylabel=w_label,
                 x_line=w_bins_avg_magnitudes,
                 y_line=avg_velocities_w,
                 yerr=velocities_w_std,
                 x_scatter=w_magnitudes,
                 y_scatter=w_velocities)

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
                 x_line: List[float],
                 y_line: List[float],
                 yerr: List[float],
                 marker: str = 's',
                 markersize: float = 3.,
                 line_color: str = 'k',
                 capsize: float = 5.,
                 linewidth: float = 1.,
                 x_scatter: List[float],
                 y_scatter: List[float],
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

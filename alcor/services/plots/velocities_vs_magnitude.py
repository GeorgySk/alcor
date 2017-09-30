from typing import (Tuple,
                    List)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sqlalchemy.orm.session import Session

from alcor.models.velocities_vs_magnitudes.bins import (Bin,
                                                        LepineCaseUBin,
                                                        LepineCaseVBin,
                                                        LepineCaseWBin)
from alcor.models.velocities_vs_magnitudes.clouds import (Cloud,
                                                          LepineCaseUCloud,
                                                          LepineCaseVCloud,
                                                          LepineCaseWCloud)
from alcor.services.data_access import fetch_all


def plot(session: Session,
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

    # TODO: implement other ways of fetching
    bins = fetch_all(Bin,
                     session=session)

    avg_bin_magnitudes = [stars_bin.avg_magnitude
                          for stars_bin in bins]
    avg_u_velocities = [stars_bin.avg_u_velocity
                        for stars_bin in bins]
    avg_velocities_v = [stars_bin.avg_v_velocity
                        for stars_bin in bins]
    avg_velocities_w = [stars_bin.avg_w_velocity
                        for stars_bin in bins]
    u_velocities_std = [stars_bin.u_velocity_std
                        for stars_bin in bins]
    velocities_v_std = [stars_bin.v_velocity_std
                        for stars_bin in bins]
    velocities_w_std = [stars_bin.w_velocity_std
                        for stars_bin in bins]

    # TODO: implement other ways of fetching
    clouds = fetch_all(Cloud,
                       session=session)

    magnitudes = [star.bolometric_magnitude
                  for star in clouds]
    u_velocities = [star.u_velocity
                    for star in clouds]
    v_velocities = [star.v_velocity
                    for star in clouds]
    w_velocities = [star.w_velocity
                    for star in clouds]

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


def plot_lepine_case(session: Session,
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

    # TODO: implement other fetching functions
    u_vs_mag_bins = fetch_all(LepineCaseUBin,
                              session=session)
    v_vs_mag_bins = fetch_all(LepineCaseVBin,
                              session=session)
    w_vs_mag_bins = fetch_all(LepineCaseWBin,
                              session=session)

    u_bins_avg_magnitudes = [stars_bin.avg_magnitude
                             for stars_bin in u_vs_mag_bins]
    avg_u_velocities = [stars_bin.avg_u_velocity
                        for stars_bin in u_vs_mag_bins]
    u_velocities_std = [stars_bin.u_velocity_std
                        for stars_bin in u_vs_mag_bins]
    v_bins_avg_magnitudes = [stars_bin.avg_magnitude
                             for stars_bin in v_vs_mag_bins]
    avg_velocities_v = [stars_bin.avg_v_velocity
                        for stars_bin in v_vs_mag_bins]
    velocities_v_std = [stars_bin.v_velocity_std
                        for stars_bin in v_vs_mag_bins]
    w_bins_avg_magnitudes = [stars_bin.avg_magnitude
                             for stars_bin in w_vs_mag_bins]
    avg_velocities_w = [stars_bin.avg_w_velocity
                        for stars_bin in w_vs_mag_bins]
    velocities_w_std = [stars_bin.w_velocity_std
                        for stars_bin in w_vs_mag_bins]

    # TODO: implement other fetching functions
    u_vs_mag_cloud = fetch_all(LepineCaseUCloud,
                               session=session)
    v_vs_mag_cloud = fetch_all(LepineCaseVCloud,
                               session=session)
    w_vs_mag_cloud = fetch_all(LepineCaseWCloud,
                               session=session)

    u_magnitudes = [star.bolometric_magnitude
                    for star in u_vs_mag_cloud]
    u_velocities = [star.u_velocity
                    for star in u_vs_mag_cloud]
    v_magnitudes = [star.bolometric_magnitude
                    for star in v_vs_mag_cloud]
    v_velocities = [star.v_velocity
                    for star in v_vs_mag_cloud]
    w_magnitudes = [star.bolometric_magnitude
                    for star in w_vs_mag_cloud]
    w_velocities = [star.w_velocity
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

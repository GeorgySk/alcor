import os
from functools import partial
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
import pandas as pd

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

    (avg_bin_magnitudes,
     avg_u_velocities,
     avg_velocities_v,
     avg_velocities_w,
     u_velocities_std,
     velocities_v_std,
     velocities_w_std) = zip(*sorted(zip(avg_bin_magnitudes,
                                         avg_u_velocities,
                                         avg_velocities_v,
                                         avg_velocities_w,
                                         u_velocities_std,
                                         velocities_v_std,
                                         velocities_w_std)))

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
                     magnitude_label: str = '$M_{bol}$',
                     u_observational_bins_path: str = 'observational_data'
                                                      '/mbol_avg_u.dat',
                     v_observational_bins_path: str = 'observational_data'
                                                      '/mbol_avg_v.dat',
                     w_observational_bins_path: str = 'observational_data'
                                                      '/mbol_avg_w.dat',
                     u_observational_clouds_path: str = 'observational_data'
                                                        '/mbol_cloud_u.dat',
                     v_observational_clouds_path: str = 'observational_data'
                                                        '/mbol_cloud_v.dat',
                     w_observational_clouds_path: str = 'observational_data'
                                                        '/mbol_cloud_w.dat',
                     ) -> None:
    subplots = {}
    figure, (subplots['u'],
             subplots['v'],
             subplots['w']) = plt.subplots(nrows=3,
                                           figsize=figure_size)

    velocities_labels = dict(u=u_label,
                             v=v_label,
                             w=w_label)
    magnitude_labels = dict(u=None,
                            v=None,
                            w=magnitude_label)

    # TODO: find a proper way to use relative paths
    base_dir = os.path.dirname(__file__)
    observational_bins_paths = dict(
            u=os.path.join(base_dir, u_observational_bins_path),
            v=os.path.join(base_dir, v_observational_bins_path),
            w=os.path.join(base_dir, w_observational_bins_path))
    observational_clouds_paths = dict(
            u=os.path.join(base_dir, u_observational_clouds_path),
            v=os.path.join(base_dir, v_observational_clouds_path),
            w=os.path.join(base_dir, w_observational_clouds_path))

    read_csv = partial(pd.read_csv,
                       delimiter=' ',
                       skipinitialspace=True,
                       header=None)

    # TODO: implement other fetching functions
    fetch = partial(fetch_all,
                    session=session)
    velocity_vs_magnitude_bins = dict(u=fetch(LepineCaseUBin),
                                      v=fetch(LepineCaseVBin),
                                      w=fetch(LepineCaseWBin))
    velocity_vs_magnitude_clouds = dict(u=fetch(LepineCaseUCloud),
                                        v=fetch(LepineCaseVCloud),
                                        w=fetch(LepineCaseWCloud))

    observational_bins_dfs = {}
    observational_clouds_dfs = {}
    bins_avg_magnitudes = {}
    magnitudes = {}
    for key in ('u', 'v', 'w'):
        observational_bins_dfs[key] = read_csv(observational_bins_paths[key])
        observational_clouds_dfs[key] = read_csv(
                observational_clouds_paths[key])
        bins_avg_magnitudes[key] = [
            stars_bin.avg_magnitude
            for stars_bin in velocity_vs_magnitude_bins[key]]
        magnitudes[key] = [star.bolometric_magnitude
                           for star in velocity_vs_magnitude_clouds[key]]

    avg_velocities = dict(
            u=[stars_bin.avg_u_velocity
               for stars_bin in velocity_vs_magnitude_bins['u']],
            v=[stars_bin.avg_v_velocity
               for stars_bin in velocity_vs_magnitude_bins['v']],
            w=[stars_bin.avg_w_velocity
               for stars_bin in velocity_vs_magnitude_bins['w']])
    velocities_std = dict(
            u=[stars_bin.u_velocity_std
               for stars_bin in velocity_vs_magnitude_bins['u']],
            v=[stars_bin.v_velocity_std
               for stars_bin in velocity_vs_magnitude_bins['v']],
            w=[stars_bin.w_velocity_std
               for stars_bin in velocity_vs_magnitude_bins['w']])
    velocities = dict(u=[star.u_velocity
                         for star in velocity_vs_magnitude_clouds['u']],
                      v=[star.v_velocity
                         for star in velocity_vs_magnitude_clouds['v']],
                      w=[star.w_velocity
                         for star in velocity_vs_magnitude_clouds['w']])

    for key in ('u', 'v', 'w'):
        (bins_avg_magnitudes[key],
         avg_velocities[key],
         velocities_std[key]) = zip(*sorted(zip(bins_avg_magnitudes[key],
                                                avg_velocities[key],
                                                velocities_std[key])))
        draw_subplot(subplot=subplots[key],
                     xlabel=magnitude_labels[key],
                     ylabel=velocities_labels[key],
                     x_line=bins_avg_magnitudes[key],
                     y_line=avg_velocities[key],
                     yerr=velocities_std[key],
                     x_scatter=magnitudes[key],
                     y_scatter=velocities[key],
                     x_line_obs=observational_bins_dfs[key][0],
                     y_line_obs=observational_bins_dfs[key][1],
                     yerr_obs=observational_bins_dfs[key][2],
                     x_scatter_obs=observational_clouds_dfs[key][0],
                     y_scatter_obs=observational_clouds_dfs[key][1])

    # Removing unnecessary x-labels for top and middle subplots
    subplots['u'].set_xticklabels([])
    subplots['v'].set_xticklabels([])

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
                 x_line_obs: pd.core.series.Series = None,
                 y_line_obs: pd.core.series.Series = None,
                 yerr_obs: pd.core.series.Series = None,
                 x_scatter_obs: pd.core.series.Series = None,
                 y_scatter_obs: pd.core.series.Series = None,
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
    subplot.errorbar(x=x_line_obs,
                     y=y_line_obs,
                     yerr=yerr_obs,
                     marker=marker,
                     markersize=markersize,
                     color='r',
                     capsize=capsize,
                     linewidth=linewidth)
    subplot.scatter(x=x_scatter,
                    y=y_scatter,
                    color=scatter_color,
                    s=scatter_point_size)
    subplot.scatter(x=x_scatter_obs,
                    y=y_scatter_obs,
                    color='r',
                    s=scatter_point_size)

    subplot.minorticks_on()
    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')
    subplot.set_aspect(ratio / subplot.get_data_ratio())

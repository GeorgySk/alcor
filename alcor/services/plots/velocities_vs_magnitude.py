from typing import (Tuple,
                    List)

from matplotlib.axes import Axes
from sqlalchemy.orm.session import Session
import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.services.data_access import (fetch_all_bins,
                                        fetch_all_u_vs_mag_clouds,
                                        fetch_all_v_vs_mag_clouds,
                                        fetch_all_w_vs_mag_clouds,
                                        fetch_all_u_vs_mag_bins,
                                        fetch_all_v_vs_mag_bins,
                                        fetch_all_w_vs_mag_bins)
from alcor.models.velocities_vs_magnitudes.clouds import Cloud


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
    bins = fetch_all_bins(session=session)

    avg_bin_magnitudes = [stars_bin.avg_magnitude
                          for stars_bin in bins]
    avg_velocities_u = [stars_bin.avg_velocity_u
                        for stars_bin in bins]
    avg_velocities_v = [stars_bin.avg_velocity_v
                        for stars_bin in bins]
    avg_velocities_w = [stars_bin.avg_velocity_w
                        for stars_bin in bins]
    velocities_u_std = [stars_bin.velocity_u_std
                        for stars_bin in bins]
    velocities_v_std = [stars_bin.velocity_v_std
                        for stars_bin in bins]
    velocities_w_std = [stars_bin.velocity_w_std
                        for stars_bin in bins]

    (avg_bin_magnitudes,
     avg_velocities_u,
     avg_velocities_v,
     avg_velocities_w,
     velocities_u_std,
     velocities_v_std,
     velocities_w_std) = zip(*sorted(zip(avg_bin_magnitudes,
                                         avg_velocities_u,
                                         avg_velocities_v,
                                         avg_velocities_w,
                                         velocities_u_std,
                                         velocities_v_std,
                                         velocities_w_std)))

    # TODO: implement other ways of fetching
    clouds = fetch_all_clouds(session=session)

    magnitudes = [star.bolometric_magnitude
                  for star in clouds]
    velocities_u = [star.velocity_u
                    for star in clouds]
    velocities_v = [star.velocity_v
                    for star in clouds]
    velocities_w = [star.velocity_w
                    for star in clouds]

    draw_subplot(subplot=subplot_u,
                 ylabel=u_label,
                 x_line=avg_bin_magnitudes,
                 y_line=avg_velocities_u,
                 yerr=velocities_u_std,
                 x_scatter=magnitudes,
                 y_scatter=velocities_u)
    draw_subplot(subplot=subplot_v,
                 ylabel=v_label,
                 x_line=avg_bin_magnitudes,
                 y_line=avg_velocities_v,
                 yerr=velocities_v_std,
                 x_scatter=magnitudes,
                 y_scatter=velocities_v)
    draw_subplot(subplot=subplot_w,
                 xlabel=magnitude_label,
                 ylabel=w_label,
                 x_line=avg_bin_magnitudes,
                 y_line=avg_velocities_w,
                 yerr=velocities_w_std,
                 x_scatter=magnitudes,
                 y_scatter=velocities_w)

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
    u_vs_mag_bins = fetch_all_u_vs_mag_bins(session=session)
    v_vs_mag_bins = fetch_all_v_vs_mag_bins(session=session)
    w_vs_mag_bins = fetch_all_w_vs_mag_bins(session=session)

    u_bins_avg_magnitudes = [stars_bin.avg_magnitude
                             for stars_bin in u_vs_mag_bins]
    avg_velocities_u = [stars_bin.avg_velocity_u
                        for stars_bin in u_vs_mag_bins]
    velocities_u_std = [stars_bin.velocity_u_std
                        for stars_bin in u_vs_mag_bins]
    v_bins_avg_magnitudes = [stars_bin.avg_magnitude
                             for stars_bin in v_vs_mag_bins]
    avg_velocities_v = [stars_bin.avg_velocity_v
                        for stars_bin in v_vs_mag_bins]
    velocities_v_std = [stars_bin.velocity_v_std
                        for stars_bin in v_vs_mag_bins]
    w_bins_avg_magnitudes = [stars_bin.avg_magnitude
                             for stars_bin in w_vs_mag_bins]
    avg_velocities_w = [stars_bin.avg_velocity_w
                        for stars_bin in w_vs_mag_bins]
    velocities_w_std = [stars_bin.velocity_w_std
                        for stars_bin in w_vs_mag_bins]

    (u_bins_avg_magnitudes,
     avg_velocities_u,
     velocities_u_std) = zip(*sorted(zip(u_bins_avg_magnitudes,
                                         avg_velocities_u,
                                         velocities_u_std)))
    (v_bins_avg_magnitudes,
     avg_velocities_v,
     velocities_v_std) = zip(*sorted(zip(v_bins_avg_magnitudes,
                                         avg_velocities_v,
                                         velocities_v_std)))
    (w_bins_avg_magnitudes,
     avg_velocities_w,
     velocities_w_std) = zip(*sorted(zip(w_bins_avg_magnitudes,
                                         avg_velocities_w,
                                         velocities_w_std)))

    # TODO: implement other fetching functions
    u_vs_mag_cloud = fetch_all_u_vs_mag_clouds(session=session)
    v_vs_mag_cloud = fetch_all_v_vs_mag_clouds(session=session)
    w_vs_mag_cloud = fetch_all_w_vs_mag_clouds(session=session)

    u_magnitudes = [star.bolometric_magnitude
                    for star in u_vs_mag_cloud]
    velocities_u = [star.velocity_u
                    for star in u_vs_mag_cloud]
    v_magnitudes = [star.bolometric_magnitude
                    for star in v_vs_mag_cloud]
    velocities_v = [star.velocity_v
                    for star in v_vs_mag_cloud]
    w_magnitudes = [star.bolometric_magnitude
                    for star in w_vs_mag_cloud]
    velocities_w = [star.velocity_w
                    for star in w_vs_mag_cloud]

    draw_subplot(subplot=subplot_u,
                 ylabel=u_label,
                 x_line=u_bins_avg_magnitudes,
                 y_line=avg_velocities_u,
                 yerr=velocities_u_std,
                 x_scatter=u_magnitudes,
                 y_scatter=velocities_u)
    draw_subplot(subplot=subplot_v,
                 ylabel=v_label,
                 x_line=v_bins_avg_magnitudes,
                 y_line=avg_velocities_v,
                 yerr=velocities_v_std,
                 x_scatter=v_magnitudes,
                 y_scatter=velocities_v)
    draw_subplot(subplot=subplot_w,
                 xlabel=magnitude_label,
                 ylabel=w_label,
                 x_line=w_bins_avg_magnitudes,
                 y_line=avg_velocities_w,
                 yerr=velocities_w_std,
                 x_scatter=w_magnitudes,
                 y_scatter=velocities_w)

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


# TODO: model name collision
def fetch_all_clouds(session: Session) -> List[Cloud]:
    query = session.query(Cloud)
    return query.all()

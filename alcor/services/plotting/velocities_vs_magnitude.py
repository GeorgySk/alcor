from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.models.velocities_vs_magnitudes.bins import (Bin,
                                                        LepineCaseUBin,
                                                        LepineCaseVBin,
                                                        LepineCaseWBin)
from alcor.models.velocities_vs_magnitudes.clouds import (Cloud,
                                                          LepineCaseUCloud,
                                                          LepineCaseVCloud,
                                                          LepineCaseWCloud)
from alcor.services.data_access.reading import fetch
from alcor.utils import get_columns

FILENAME = 'velocities_vs_magnitude.ps'

# TODO: figure out how all these dimensions work
FIGURE_SIZE = (10, 12)
DESIRED_DIMENSIONS_RATIO = 7 / 13

U_LABEL = '$U_{LSR}(km/s)$'
V_LABEL = '$V_{LSR}(km/s)$'
W_LABEL = '$W_{LSR}(km/s)$'
MAGNITUDE_LABEL = '$M_{bol}$'

MAGNITUDE_LIMITS = [6, 19]
VELOCITIES_LIMITS = [-150, 150]

LINE_COLOR = 'k'
CLOUD_COLOR = 'gray'
MARKER = 's'
MARKERSIZE = 3
CAP_SIZE = 5
LINEWIDTH = 1
CLOUD_POINT_SIZE = 1

PLOT_WIDTH = 500
PLOT_HEIGHT = 250

CSV_DELIMITER = ' '


def plot(session: Session) -> None:

    # TODO: fetch only last by time(ok?) bins
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
     velocities_w_std) = (_ for _ in zip(*sorted(zip(avg_bin_magnitudes,
                                                     avg_velocities_u,
                                                     avg_velocities_v,
                                                     avg_velocities_w,
                                                     velocities_u_std,
                                                     velocities_v_std,
                                                     velocities_w_std))))

    # TODO: do I need to use sharex or sharey attrs?
    figure, (subplot_u,
             subplot_v,
             subplot_w) = plt.subplots(nrows=3,
                                       figsize=FIGURE_SIZE)

    # TODO: find the way to apply limits once for all subplots
    subplot_u.set(ylabel=U_LABEL,
                  xlim=MAGNITUDE_LIMITS,
                  ylim=VELOCITIES_LIMITS)
    subplot_v.set(ylabel=V_LABEL,
                  xlim=MAGNITUDE_LIMITS,
                  ylim=VELOCITIES_LIMITS)
    subplot_w.set(xlabel=MAGNITUDE_LABEL,
                  ylabel=W_LABEL,
                  xlim=MAGNITUDE_LIMITS,
                  ylim=VELOCITIES_LIMITS)

    subplot_u.errorbar(x=avg_bin_magnitudes,
                       y=avg_velocities_u,
                       yerr=velocities_u_std,
                       marker=MARKER,
                       markersize=MARKERSIZE,
                       color=LINE_COLOR,
                       capsize=CAP_SIZE,
                       linewidth=LINEWIDTH)
    subplot_v.errorbar(x=avg_bin_magnitudes,
                       y=avg_velocities_v,
                       yerr=velocities_v_std,
                       marker=MARKER,
                       markersize=MARKERSIZE,
                       color=LINE_COLOR,
                       capsize=CAP_SIZE,
                       linewidth=LINEWIDTH)
    subplot_w.errorbar(x=avg_bin_magnitudes,
                       y=avg_velocities_w,
                       yerr=velocities_w_std,
                       marker=MARKER,
                       markersize=MARKERSIZE,
                       color=LINE_COLOR,
                       capsize=CAP_SIZE,
                       linewidth=LINEWIDTH)

    # TODO: fetch only last by time(ok?) clouds
    clouds = fetch_all_clouds(session=session)

    magnitudes = [star.bolometric_magnitude
                  for star in clouds]
    velocities_u = [star.velocity_u
                    for star in clouds]
    velocities_v = [star.velocity_v
                    for star in clouds]
    velocities_w = [star.velocity_w
                    for star in clouds]

    subplot_u.scatter(x=magnitudes,
                      y=velocities_u,
                      color=CLOUD_COLOR,
                      s=CLOUD_POINT_SIZE)
    subplot_v.scatter(x=magnitudes,
                      y=velocities_v,
                      color=CLOUD_COLOR,
                      s=CLOUD_POINT_SIZE)
    subplot_w.scatter(x=magnitudes,
                      y=velocities_w,
                      color=CLOUD_COLOR,
                      s=CLOUD_POINT_SIZE)

    # TODO: why does this apply minorticks only to the last subplot?
    plt.minorticks_on()

    subplot_u.xaxis.set_ticks_position('both')
    subplot_u.yaxis.set_ticks_position('both')
    subplot_v.xaxis.set_ticks_position('both')
    subplot_v.yaxis.set_ticks_position('both')
    subplot_w.xaxis.set_ticks_position('both')
    subplot_w.yaxis.set_ticks_position('both')

    subplot_u.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot_u.get_data_ratio())
    subplot_v.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot_v.get_data_ratio())
    subplot_w.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot_w.get_data_ratio())

    # TODO: delete overlapping y-labels
    # TODO: delete unnecessary x-labels for top and middle subplots
    figure.subplots_adjust(hspace=0)

    # FIXME: cloud and bins are not correlated!
    plt.savefig(FILENAME)


def plot_lepine_case(session: Session) -> None:

    # TODO: fetch only last by time(ok?) bins
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
     velocities_u_std) = (_ for _ in zip(*sorted(zip(u_bins_avg_magnitudes,
                                                     avg_velocities_u,
                                                     velocities_u_std))))
    (v_bins_avg_magnitudes,
     avg_velocities_v,
     velocities_v_std) = (_ for _ in zip(*sorted(zip(v_bins_avg_magnitudes,
                                                     avg_velocities_v,
                                                     velocities_v_std))))
    (w_bins_avg_magnitudes,
     avg_velocities_w,
     velocities_w_std) = (_ for _ in zip(*sorted(zip(w_bins_avg_magnitudes,
                                                     avg_velocities_w,
                                                     velocities_w_std))))

    # TODO: do I need to use sharex or sharey attrs?
    figure, (subplot_u,
             subplot_v,
             subplot_w) = plt.subplots(nrows=3,
                                       figsize=FIGURE_SIZE)

    # TODO: find the way to apply limits once for all subplots
    subplot_u.set(ylabel=U_LABEL,
                  xlim=MAGNITUDE_LIMITS,
                  ylim=VELOCITIES_LIMITS)
    subplot_v.set(ylabel=V_LABEL,
                  xlim=MAGNITUDE_LIMITS,
                  ylim=VELOCITIES_LIMITS)
    subplot_w.set(xlabel=MAGNITUDE_LABEL,
                  ylabel=W_LABEL,
                  xlim=MAGNITUDE_LIMITS,
                  ylim=VELOCITIES_LIMITS)

    subplot_u.errorbar(x=u_bins_avg_magnitudes,
                       y=avg_velocities_u,
                       yerr=velocities_u_std,
                       marker=MARKER,
                       markersize=MARKERSIZE,
                       color=LINE_COLOR,
                       capsize=CAP_SIZE,
                       linewidth=LINEWIDTH)
    subplot_v.errorbar(x=v_bins_avg_magnitudes,
                       y=avg_velocities_v,
                       yerr=velocities_v_std,
                       marker=MARKER,
                       markersize=MARKERSIZE,
                       color=LINE_COLOR,
                       capsize=CAP_SIZE,
                       linewidth=LINEWIDTH)
    subplot_w.errorbar(x=w_bins_avg_magnitudes,
                       y=avg_velocities_w,
                       yerr=velocities_w_std,
                       marker=MARKER,
                       markersize=MARKERSIZE,
                       color=LINE_COLOR,
                       capsize=CAP_SIZE,
                       linewidth=LINEWIDTH)

    # TODO: fetch only last by time(ok?) clouds
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

    subplot_u.scatter(x=u_magnitudes,
                      y=velocities_u,
                      color=CLOUD_COLOR,
                      s=CLOUD_POINT_SIZE)
    subplot_v.scatter(x=v_magnitudes,
                      y=velocities_v,
                      color=CLOUD_COLOR,
                      s=CLOUD_POINT_SIZE)
    subplot_w.scatter(x=w_magnitudes,
                      y=velocities_w,
                      color=CLOUD_COLOR,
                      s=CLOUD_POINT_SIZE)

    # TODO: why does this apply minorticks only to the last subplot?
    plt.minorticks_on()

    subplot_u.xaxis.set_ticks_position('both')
    subplot_u.yaxis.set_ticks_position('both')
    subplot_v.xaxis.set_ticks_position('both')
    subplot_v.yaxis.set_ticks_position('both')
    subplot_w.xaxis.set_ticks_position('both')
    subplot_w.yaxis.set_ticks_position('both')

    subplot_u.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot_u.get_data_ratio())
    subplot_v.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot_v.get_data_ratio())
    subplot_w.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot_w.get_data_ratio())

    # TODO: delete overlapping y-labels
    # TODO: delete unnecessary x-labels for top and middle subplots
    figure.subplots_adjust(hspace=0)

    # FIXME: cloud and bins are not correlated!
    plt.savefig(FILENAME)


def fetch_all_bins(*,
                   session: Session):
    query = (Bin.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [Bin(**record)
            for record in records]


def fetch_all_clouds(*,
                     session: Session):
    query = (Cloud.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [Cloud(**record)
            for record in records]


def fetch_all_u_vs_mag_bins(*,
                            session: Session):
    query = (LepineCaseUBin.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [LepineCaseUBin(**record)
            for record in records]


def fetch_all_v_vs_mag_bins(*,
                            session: Session):
    query = (LepineCaseVBin.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [LepineCaseVBin(**record)
            for record in records]


def fetch_all_w_vs_mag_bins(*,
                            session: Session):
    query = (LepineCaseWBin.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [LepineCaseWBin(**record)
            for record in records]


def fetch_all_u_vs_mag_clouds(*,
                              session: Session):
    query = (LepineCaseUCloud.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [LepineCaseUCloud(**record)
            for record in records]


def fetch_all_v_vs_mag_clouds(*,
                              session: Session):
    query = (LepineCaseVCloud.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [LepineCaseVCloud(**record)
            for record in records]


def fetch_all_w_vs_mag_clouds(*,
                              session: Session):
    query = (LepineCaseWCloud.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [LepineCaseWCloud(**record)
            for record in records]

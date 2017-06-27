import csv
from typing import List

from bokeh.io import save
from bokeh.layouts import gridplot
from bokeh.plotting import (figure,
                            output_file)
from bokeh.plotting.figure import Figure
from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.models.velocities_vs_magnitudes.bins import Bin
from alcor.models.velocities_vs_magnitudes.clouds import Cloud
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

    avg_bin_magnitudes = [_.avg_magnitude
                          for _ in bins]
    avg_velocities_u = [_.avg_velocity_u
                        for _ in bins]
    avg_velocities_v = [_.avg_velocity_v
                        for _ in bins]
    avg_velocities_w = [_.avg_velocity_w
                        for _ in bins]
    velocities_u_std = [_.velocity_u_std
                        for _ in bins]
    velocities_v_std = [_.velocity_v_std
                        for _ in bins]
    velocities_w_std = [_.velocity_w_std
                        for _ in bins]

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

    magnitudes = [_.bolometric_magnitude
                  for _ in clouds]
    velocities_u = [_.velocity_u
                    for _ in clouds]
    velocities_v = [_.velocity_v
                    for _ in clouds]
    velocities_w = [_.velocity_w
                    for _ in clouds]

    (magnitudes,
     velocities_u,
     velocities_v,
     velocities_w) = (_ for _ in zip(*sorted(zip(magnitudes,
                                                 velocities_u,
                                                 velocities_v,
                                                 velocities_w))))

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


def plot_lepine_case(session: Session):
    output_file("velocities_vs_magnitude.html")

    top_plot = figure()
    middle_plot = figure()
    bottom_plot = figure()

    with open('u_vs_mag_bins.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (average_bin_magnitude_u,
         average_velocity_u,
         velocity_u_std) = get_columns(reader)
    top_plot.line(average_bin_magnitude_u,
                  average_velocity_u)
    top_plot.square(average_bin_magnitude_u,
                    average_velocity_u)
    add_errorbars(fig=top_plot,
                  x=average_bin_magnitude_u,
                  y=average_velocity_u,
                  y_err=velocity_u_std)

    with open('v_vs_mag_bins.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (average_bin_magnitude_v,
         average_velocity_v,
         velocity_v_std) = get_columns(reader)
    middle_plot.line(average_bin_magnitude_v,
                     average_velocity_v)
    middle_plot.square(average_bin_magnitude_v,
                       average_velocity_v)
    add_errorbars(fig=middle_plot,
                  x=average_bin_magnitude_v,
                  y=average_velocity_v,
                  y_err=velocity_v_std)

    with open('w_vs_mag_bins.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (average_bin_magnitude_w,
         average_velocity_w,
         velocity_w_std) = get_columns(reader)
    bottom_plot.line(average_bin_magnitude_w,
                     average_velocity_w)
    bottom_plot.square(average_bin_magnitude_w,
                       average_velocity_w)
    add_errorbars(fig=bottom_plot,
                  x=average_bin_magnitude_w,
                  y=average_velocity_w,
                  y_err=velocity_w_std)

    with open('u_vs_mag_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (bolometric_magnitude,
         velocity_u) = get_columns(reader)
    top_plot.circle(bolometric_magnitude,
                    velocity_u,
                    size=1)

    with open('v_vs_mag_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (bolometric_magnitude,
         velocity_v) = get_columns(reader)
    middle_plot.circle(bolometric_magnitude,
                       velocity_v,
                       size=1)

    with open('w_vs_mag_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (bolometric_magnitude,
         velocity_w) = get_columns(reader)
    bottom_plot.circle(bolometric_magnitude,
                       velocity_w,
                       size=1)

    main_plot = gridplot(children=[top_plot,
                                   middle_plot,
                                   bottom_plot],
                         ncols=1,
                         plot_width=PLOT_WIDTH,
                         plot_height=PLOT_HEIGHT,
                         merge_tools=True,
                         toolbar_location='right')
    save(main_plot)


def add_errorbars(fig: Figure,
                  x: List[float],
                  y: List[float],
                  y_err: List[float]) -> None:
    multiline_x = []
    multiline_y = []
    for (x_coordinate, y_coordinate, error) in zip(x, y, y_err):
        multiline_x.append((x_coordinate,
                            x_coordinate))
        multiline_y.append((y_coordinate - error,
                            y_coordinate + error))

    fig.multi_line(multiline_x,
                   multiline_y)


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

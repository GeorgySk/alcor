import csv
from typing import List

from bokeh.io import save
from bokeh.layouts import gridplot
from bokeh.plotting import (figure,
                            output_file)
from bokeh.plotting.figure import Figure

from alcor.utils import get_columns


PLOT_WIDTH = 500
PLOT_HEIGHT = 250

CSV_DELIMITER = ' '


def plot() -> None:
    output_file("velocities_vs_magnitude.html")

    # TODO: maybe I should use Plot here
    top_plot = figure()
    middle_plot = figure()
    bottom_plot = figure()

    with open('uvw_vs_mag_bins.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (average_bin_magnitude,
         average_velocity_u,
         average_velocity_v,
         average_velocity_w,
         velocity_u_std,
         velocity_v_std,
         velocity_w_std) = get_columns(reader)

    top_plot.line(average_bin_magnitude,
                  average_velocity_u)
    top_plot.square(average_bin_magnitude,
                    average_velocity_u)

    middle_plot.line(average_bin_magnitude,
                     average_velocity_v)
    middle_plot.square(average_bin_magnitude,
                       average_velocity_v)

    bottom_plot.line(average_bin_magnitude,
                     average_velocity_w)
    bottom_plot.square(average_bin_magnitude,
                       average_velocity_w)

    add_errorbars(fig=top_plot,
                  x=average_bin_magnitude,
                  y=average_velocity_u,
                  y_err=velocity_u_std)
    add_errorbars(fig=middle_plot,
                  x=average_bin_magnitude,
                  y=average_velocity_v,
                  y_err=velocity_v_std)
    add_errorbars(fig=bottom_plot,
                  x=average_bin_magnitude,
                  y=average_velocity_w,
                  y_err=velocity_w_std)

    with open('uvw_vs_mag_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (bolometric_magnitude,
         velocity_u,
         velocity_v,
         velocity_w) = get_columns(reader)

    top_plot.circle(bolometric_magnitude,
                    velocity_u,
                    size=1)
    middle_plot.circle(bolometric_magnitude,
                       velocity_v,
                       size=1)
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


def plot_lepine_case():
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

    with open('v_vs_mag_bins.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (average_bin_magnitude_v,
         average_velocity_v,
         velocity_v_std) = get_columns(reader)

    with open('w_vs_mag_bins.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (average_bin_magnitude_w,
         average_velocity_w,
         velocity_w_std) = get_columns(reader)

    top_plot.line(average_bin_magnitude_u,
                  average_velocity_u)
    top_plot.square(average_bin_magnitude_u,
                    average_velocity_u)

    middle_plot.line(average_bin_magnitude_v,
                     average_velocity_v)
    middle_plot.square(average_bin_magnitude_v,
                       average_velocity_v)

    bottom_plot.line(average_bin_magnitude_w,
                     average_velocity_w)
    bottom_plot.square(average_bin_magnitude_w,
                       average_velocity_w)

    add_errorbars(fig=top_plot,
                  x=average_bin_magnitude_u,
                  y=average_velocity_u,
                  y_err=velocity_u_std)
    add_errorbars(fig=middle_plot,
                  x=average_bin_magnitude_v,
                  y=average_velocity_v,
                  y_err=velocity_v_std)
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

    with open('v_vs_mag_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (bolometric_magnitude,
         velocity_v) = get_columns(reader)

    with open('w_vs_mag_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (bolometric_magnitude,
         velocity_w) = get_columns(reader)

    top_plot.circle(bolometric_magnitude,
                    velocity_u,
                    size=1)
    middle_plot.circle(bolometric_magnitude,
                       velocity_v,
                       size=1)
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

import csv

from bokeh.layouts import column
from bokeh.plotting import (figure,
                            output_file,
                            show)


def plot() -> None:
    output_file("velocities_vs_magnitude.html")

    with open('uvw_vs_mag_bins.csv', 'r') as data_file:
        reader = csv.reader(data_file,
                            delimiter=' ')

        header_row = next(reader)

        average_bin_magnitude = []
        average_velocity_u = []
        average_velocity_v = []
        average_velocity_w = []
        velocity_u_std = []
        velocity_v_std = []
        velocity_w_std = []

        for row in reader:
            average_bin_magnitude.append(float(row[0]))
            average_velocity_u.append(float(row[1]))
            average_velocity_v.append(float(row[2]))
            average_velocity_w.append(float(row[3]))
            velocity_u_std.append(float(row[4]))
            velocity_v_std.append(float(row[5]))
            velocity_w_std.append(float(row[6]))

    top_plot = figure(width=500, height=250)
    top_plot.line(average_bin_magnitude,
                  average_velocity_u)
    top_plot.square(average_bin_magnitude,
                    average_velocity_u)

    middle_plot = figure(width=500, height=250)
    middle_plot.line(average_bin_magnitude,
                     average_velocity_v)
    middle_plot.square(average_bin_magnitude,
                       average_velocity_v)

    bottom_plot = figure(width=500, height=250)
    bottom_plot.line(average_bin_magnitude,
                     average_velocity_w)
    bottom_plot.square(average_bin_magnitude,
                       average_velocity_w)

    multiline_x = []
    multiline_y_u = []
    multiline_y_v = []
    multiline_y_w = []

    for (magnitude,
         velocity_u, velocity_v, velocity_w,
         std_u, std_v, std_w) in zip(average_bin_magnitude,
                                     average_velocity_u,
                                     average_velocity_v,
                                     average_velocity_w,
                                     velocity_u_std,
                                     velocity_v_std,
                                     velocity_w_std):
        multiline_x.append((magnitude,
                            magnitude))
        multiline_y_u.append((velocity_u - std_u,
                              velocity_u + std_u))
        multiline_y_v.append((velocity_v - std_v,
                              velocity_v + std_v))
        multiline_y_w.append((velocity_w - std_w,
                              velocity_w + std_w))

    top_plot.multi_line(multiline_x,
                        multiline_y_u)
    middle_plot.multi_line(multiline_x,
                           multiline_y_v)
    bottom_plot.multi_line(multiline_x,
                           multiline_y_w)

    with open('uvw_vs_mag_cloud.csv', 'r') as data_file:
        reader = csv.reader(data_file,
                            delimiter=' ')

        header_row = next(reader)

        bolometric_magnitude = []
        velocity_u = []
        velocity_v = []
        velocity_w = []

        for row in reader:
            bolometric_magnitude.append(float(row[0]))
            velocity_u.append(float(row[1]))
            velocity_v.append(float(row[2]))
            velocity_w.append(float(row[3]))

    top_plot.circle(bolometric_magnitude,
                    velocity_u,
                    size=1)
    middle_plot.circle(bolometric_magnitude,
                       velocity_v,
                       size=1)
    bottom_plot.circle(bolometric_magnitude,
                       velocity_w,
                       size=1)

    show(column(top_plot,
                middle_plot,
                bottom_plot))


def plot_lepine_case():
    pass

import csv

from bokeh.layouts import column
from bokeh.plotting import (figure,
                            output_file,
                            show)


def plot() -> None:
    output_file("velocity_clouds.html")

    with open('uvw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')

        header_row = next(reader)

        velocity_u = []
        velocity_v = []
        velocity_w = []

        for row in reader:
            velocity_u.append(float(row[0]))
            velocity_v.append(float(row[1]))
            velocity_w.append(float(row[2]))

    top_plot = figure(width=250, height=250)
    top_plot.circle(velocity_u,
                    velocity_w,
                    size=1)

    middle_plot = figure(width=250, height=250)
    middle_plot.circle(velocity_u,
                       velocity_v,
                       size=1)

    bottom_plot = figure(width=250, height=250)
    bottom_plot.circle(velocity_w,
                       velocity_v,
                       size=1)

    show(column(top_plot,
                middle_plot,
                bottom_plot))


def plot_lepine_case():
    output_file("velocity_clouds.html")

    with open('uw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')
        header_row = next(reader)
        velocity_u = []
        velocity_w = []
        for row in reader:
            velocity_u.append(float(row[0]))
            velocity_w.append(float(row[1]))
    top_plot = figure(width=250, height=250)
    top_plot.circle(velocity_u,
                    velocity_w,
                    size=1)

    with open('uv_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')
        header_row = next(reader)
        velocity_u = []
        velocity_v = []
        for row in reader:
            velocity_u.append(float(row[0]))
            velocity_v.append(float(row[1]))
    middle_plot = figure(width=250, height=250)
    middle_plot.circle(velocity_u,
                       velocity_v,
                       size=1)

    with open('vw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')
        header_row = next(reader)
        velocity_v = []
        velocity_w = []
        for row in reader:
            velocity_v.append(float(row[0]))
            velocity_w.append(float(row[1]))
    bottom_plot = figure(width=250, height=250)
    bottom_plot.circle(velocity_v,
                       velocity_w,
                       size=1)

    show(column(top_plot,
                middle_plot,
                bottom_plot))

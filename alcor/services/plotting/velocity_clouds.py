import csv

from bokeh.layouts import column
from bokeh.models.glyphs import Ellipse
from bokeh.plotting import (figure,
                            output_file,
                            show)

# Kinematic properties of the thin disk taken from the paper of
# N.Rowell and N.C.Hambly (mean motions are relative to the Sun):
# "White dwarfs in the SuperCOSMOS Sky Survey: the thin disc,
# thick disc and spheroid luminosity functions"
# Mon. Not. R. Astron. Soc. 417, 93–113 (2011)
# doi:10.1111/j.1365-2966.2011.18976.x
AVERAGE_POPULATION_VELOCITY_U = -8.62
AVERAGE_POPULATION_VELOCITY_V = -20.04
AVERAGE_POPULATION_VELOCITY_W = -7.1
STD_POPULATION_U = 32.4
STD_POPULATION_V = 23.
STD_POPULATION_W = 18.1


def plot() -> None:
    output_file("velocity_clouds.html")

    with open('uvw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')

        header_row = next(reader)

        rows = (map(float, row) for row in reader)
        (velocity_u, velocity_v, velocity_w) = zip(*rows)

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

    top_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                          y=AVERAGE_POPULATION_VELOCITY_W,
                          width=STD_POPULATION_U * 2,
                          height=STD_POPULATION_W * 2,
                          fill_color='white',
                          fill_alpha=0.,
                          line_dash='dashed')
    top_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                                 y=AVERAGE_POPULATION_VELOCITY_W,
                                 width=STD_POPULATION_U * 4,
                                 height=STD_POPULATION_W * 4,
                                 fill_color='white',
                                 fill_alpha=0.)
    top_plot.add_glyph(top_ellipse)
    top_plot.add_glyph(top_double_ellipse)

    middle_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                             y=AVERAGE_POPULATION_VELOCITY_V,
                             width=STD_POPULATION_U * 2,
                             height=STD_POPULATION_V * 2,
                             fill_color='white',
                             fill_alpha=0.,
                             line_dash='dashed')
    middle_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                                    y=AVERAGE_POPULATION_VELOCITY_V,
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_V * 4,
                                    fill_color='white',
                                    fill_alpha=0.)
    middle_plot.add_glyph(middle_ellipse)
    middle_plot.add_glyph(middle_double_ellipse)

    bottom_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
                             y=AVERAGE_POPULATION_VELOCITY_W,
                             width=STD_POPULATION_U * 2,
                             height=STD_POPULATION_V * 2,
                             fill_color='white',
                             fill_alpha=0.,
                             line_dash='dashed')
    bottom_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
                                    y=AVERAGE_POPULATION_VELOCITY_W,
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_V * 4,
                                    fill_color='white',
                                    fill_alpha=0.)
    bottom_plot.add_glyph(bottom_ellipse)
    bottom_plot.add_glyph(bottom_double_ellipse)

    show(column(top_plot,
                middle_plot,
                bottom_plot))


def plot_lepine_case():
    output_file("velocity_clouds.html")

    with open('uw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')
        header_row = next(reader)
        rows = (map(float, row) for row in reader)
        (velocity_u, velocity_w) = zip(*rows)
    top_plot = figure(width=250, height=250)
    top_plot.circle(velocity_u,
                    velocity_w,
                    size=1)

    with open('uv_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')
        header_row = next(reader)
        rows = (map(float, row) for row in reader)
        (velocity_u, velocity_v) = zip(*rows)
    middle_plot = figure(width=250, height=250)
    middle_plot.circle(velocity_u,
                       velocity_v,
                       size=1)

    with open('vw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')
        header_row = next(reader)
        rows = (map(float, row) for row in reader)
        (velocity_v, velocity_w) = zip(*rows)
    bottom_plot = figure(width=250, height=250)
    bottom_plot.circle(velocity_v,
                       velocity_w,
                       size=1)

    top_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                          y=AVERAGE_POPULATION_VELOCITY_W,
                          width=STD_POPULATION_U * 2,
                          height=STD_POPULATION_W * 2,
                          fill_color='white',
                          fill_alpha=0.,
                          line_dash='dashed')
    top_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                                 y=AVERAGE_POPULATION_VELOCITY_W,
                                 width=STD_POPULATION_U * 4,
                                 height=STD_POPULATION_W * 4,
                                 fill_color='white',
                                 fill_alpha=0.)
    top_plot.add_glyph(top_ellipse)
    top_plot.add_glyph(top_double_ellipse)

    middle_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                             y=AVERAGE_POPULATION_VELOCITY_V,
                             width=STD_POPULATION_U * 2,
                             height=STD_POPULATION_V * 2,
                             fill_color='white',
                             fill_alpha=0.)
    middle_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
                                    y=AVERAGE_POPULATION_VELOCITY_V,
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_V * 4,
                                    fill_color='white',
                                    fill_alpha=0.)
    middle_plot.add_glyph(middle_ellipse)
    middle_plot.add_glyph(middle_double_ellipse)

    bottom_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
                             y=AVERAGE_POPULATION_VELOCITY_W,
                             width=STD_POPULATION_U * 2,
                             height=STD_POPULATION_V * 2,
                             fill_color='white',
                             fill_alpha=0.)
    bottom_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
                                    y=AVERAGE_POPULATION_VELOCITY_W,
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_V * 4,
                                    fill_color='white',
                                    fill_alpha=0.)
    bottom_plot.add_glyph(bottom_ellipse)
    bottom_plot.add_glyph(bottom_double_ellipse)

    show(column(top_plot,
                middle_plot,
                bottom_plot))

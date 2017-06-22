import csv

from bokeh.io import save
from bokeh.layouts import (column,
                           gridplot)
from bokeh.models.glyphs import Ellipse
from bokeh.plotting import (figure,
                            output_file,
                            show)
from cassandra.cluster import Session

# Kinematic properties of the thin disk taken from the paper of
# N.Rowell and N.C.Hambly (mean motions are relative to the Sun):
# "White dwarfs in the SuperCOSMOS Sky Survey: the thin disc,
# thick disc and spheroid luminosity functions"
# Mon. Not. R. Astron. Soc. 417, 93–113 (2011)
# doi:10.1111/j.1365-2966.2011.18976.x
from alcor.utils import get_columns

AVERAGE_POPULATION_VELOCITY_U = -8.62
AVERAGE_POPULATION_VELOCITY_V = -20.04
AVERAGE_POPULATION_VELOCITY_W = -7.1
STD_POPULATION_U = 32.4
STD_POPULATION_V = 23.
STD_POPULATION_W = 18.1

CSV_DELIMITER = ' '

PLOT_WIDTH = 250
PLOT_HEIGHT = 250


def plot(session: Session) -> None:
    output_file("velocity_clouds.html")

    top_plot = figure()
    middle_plot = figure()
    bottom_plot = figure()

    top_plot.xaxis.axis_label = 'U'
    top_plot.yaxis.axis_label = 'W'

    middle_plot.xaxis.axis_label = 'U'
    middle_plot.yaxis.axis_label = 'V'

    bottom_plot.xaxis.axis_label = 'W'
    bottom_plot.yaxis.axis_label = 'V'

    with open('uvw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (velocity_u, velocity_v, velocity_w) = get_columns(reader)
    top_plot.circle(velocity_u,
                    velocity_w,
                    size=1)
    middle_plot.circle(velocity_u,
                       velocity_v,
                       size=1)
    bottom_plot.circle(velocity_w,
                       velocity_v,
                       size=1)

    # top_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
    #                       y=AVERAGE_POPULATION_VELOCITY_W,
    #                       width=STD_POPULATION_U * 2,
    #                       height=STD_POPULATION_W * 2,
    #                       fill_color='white',
    #                       fill_alpha=0.,
    #                       line_dash='dashed')
    # top_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
    #                              y=AVERAGE_POPULATION_VELOCITY_W,
    #                              width=STD_POPULATION_U * 4,
    #                              height=STD_POPULATION_W * 4,
    #                              fill_color='white',
    #                              fill_alpha=0.)
    # top_plot.add_glyph(top_ellipse)
    # top_plot.add_glyph(top_double_ellipse)
    #
    # middle_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
    #                          y=AVERAGE_POPULATION_VELOCITY_V,
    #                          width=STD_POPULATION_U * 2,
    #                          height=STD_POPULATION_V * 2,
    #                          fill_color='white',
    #                          fill_alpha=0.,
    #                          line_dash='dashed')
    # middle_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_U,
    #                                 y=AVERAGE_POPULATION_VELOCITY_V,
    #                                 width=STD_POPULATION_U * 4,
    #                                 height=STD_POPULATION_V * 4,
    #                                 fill_color='white',
    #                                 fill_alpha=0.)
    # middle_plot.add_glyph(middle_ellipse)
    # middle_plot.add_glyph(middle_double_ellipse)
    #
    # bottom_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
    #                          y=AVERAGE_POPULATION_VELOCITY_W,
    #                          width=STD_POPULATION_U * 2,
    #                          height=STD_POPULATION_V * 2,
    #                          fill_color='white',
    #                          fill_alpha=0.,
    #                          line_dash='dashed')
    # bottom_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
    #                                 y=AVERAGE_POPULATION_VELOCITY_W,
    #                                 width=STD_POPULATION_U * 4,
    #                                 height=STD_POPULATION_V * 4,
    #                                 fill_color='white',
    #                                 fill_alpha=0.)
    # bottom_plot.add_glyph(bottom_ellipse)
    # bottom_plot.add_glyph(bottom_double_ellipse)

    main_plot = gridplot(children=[top_plot,
                                   middle_plot,
                                   bottom_plot],
                         ncols=1,
                         plot_width=PLOT_WIDTH,
                         plot_height=PLOT_HEIGHT,
                         merge_tools=True,
                         toolbar_location='right')
    save(main_plot)


def plot_lepine_case(session: Session):
    output_file("velocity_clouds.html")

    top_plot = figure()
    middle_plot = figure()
    bottom_plot = figure()

    with open('uw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (velocity_u, velocity_w) = get_columns(reader)
    top_plot.circle(velocity_u,
                    velocity_w,
                    size=1)

    with open('uv_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (velocity_u, velocity_v) = get_columns(reader)
    middle_plot.circle(velocity_u,
                       velocity_v,
                       size=1)

    with open('vw_cloud.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        header_row = next(reader)
        (velocity_v, velocity_w) = get_columns(reader)
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
                             line_dash='dashed',
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
                             line_dash='dashed',
                             fill_alpha=0.)
    bottom_double_ellipse = Ellipse(x=AVERAGE_POPULATION_VELOCITY_V,
                                    y=AVERAGE_POPULATION_VELOCITY_W,
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_V * 4,
                                    fill_color='white',
                                    fill_alpha=0.)
    bottom_plot.add_glyph(bottom_ellipse)
    bottom_plot.add_glyph(bottom_double_ellipse)

    main_plot = gridplot(children=[top_plot,
                                   middle_plot,
                                   bottom_plot],
                         ncols=1,
                         plot_width=PLOT_WIDTH,
                         plot_height=PLOT_HEIGHT,
                         merge_tools=True,
                         toolbar_location='right')
    save(main_plot)

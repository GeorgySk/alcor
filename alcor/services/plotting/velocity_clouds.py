import csv

from bokeh.io import save
from bokeh.layouts import (column,
                           gridplot)
from bokeh.models.glyphs import Ellipse
from bokeh.plotting import (figure,
                            output_file,
                            show)
from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from alcor.utils import get_columns
from alcor.models.velocities.clouds import Cloud
from alcor.services.data_access.reading import fetch


FILENAME = 'velocity_clouds.ps'

# TODO: figure out how to work with sizes
FIGURE_SIZE = (8, 12)
DESIRED_DIMENSIONS_RATIO = 10 / 13

SUBPLOTS_SPACING = 0.25

U_LABEL = '$U(km/s)$'
V_LABEL = '$V(km/s)$'
W_LABEL = '$W(km/s)$'

U_LIMITS = [-150, 150]
V_LIMITS = [-150, 150]
W_LIMITS = [-150, 150]

CLOUD_COLOR = 'k'
POINT_SIZE = 0.5

ELLIPSE_COLOR = 'b'

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
STD_POPULATION_V = 23
STD_POPULATION_W = 18.1

CSV_DELIMITER = ' '

PLOT_WIDTH = 250
PLOT_HEIGHT = 250


def plot(session: Session) -> None:

    # TODO: Implement getting last points by time(ok?)
    cloud_points = fetch_all_cloud_points(session=session)

    velocities_u = [_.velocity_u
                    for _ in cloud_points]
    velocities_v = [_.velocity_v
                    for _ in cloud_points]
    velocities_w = [_.velocity_w
                    for _ in cloud_points]

    (velocities_u,
     velocities_v,
     velocities_w) = (_ for _ in zip(*sorted(zip(velocities_u,
                                                 velocities_v,
                                                 velocities_w))))

    # TODO: do I need to use sharex or sharey attrs?
    figure, (subplot_uv,
             subplot_uw,
             subplot_vw) = plt.subplots(nrows=3,
                                        figsize=FIGURE_SIZE)

    # TODO: find the way to apply limits once for all subplots
    subplot_uv.set(xlabel=U_LABEL,
                   ylabel=V_LABEL,
                   xlim=U_LIMITS,
                   ylim=V_LIMITS)
    subplot_uw.set(xlabel=U_LABEL,
                   ylabel=W_LABEL,
                   xlim=U_LIMITS,
                   ylim=W_LIMITS)
    subplot_vw.set(xlabel=V_LABEL,
                   ylabel=W_LABEL,
                   xlim=V_LIMITS,
                   ylim=W_LIMITS)

    subplot_uv.scatter(x=velocities_u,
                       y=velocities_v,
                       color=CLOUD_COLOR,
                       s=POINT_SIZE)
    subplot_uw.scatter(x=velocities_u,
                       y=velocities_w,
                       color=CLOUD_COLOR,
                       s=POINT_SIZE)
    subplot_vw.scatter(x=velocities_v,
                       y=velocities_w,
                       color=CLOUD_COLOR,
                       s=POINT_SIZE)

    uv_std_ellipse = Ellipse(xy=(AVERAGE_POPULATION_VELOCITY_U,
                                 AVERAGE_POPULATION_VELOCITY_V),
                             width=STD_POPULATION_U * 2,
                             height=STD_POPULATION_V * 2,
                             fill=False,
                             edgecolor=ELLIPSE_COLOR,
                             linestyle='dashed')
    uw_std_ellipse = Ellipse(xy=(AVERAGE_POPULATION_VELOCITY_U,
                                 AVERAGE_POPULATION_VELOCITY_W),
                             width=STD_POPULATION_U * 2,
                             height=STD_POPULATION_W * 2,
                             fill=False,
                             edgecolor=ELLIPSE_COLOR,
                             linestyle='dashed')
    vw_std_ellipse = Ellipse(xy=(AVERAGE_POPULATION_VELOCITY_V,
                                 AVERAGE_POPULATION_VELOCITY_W),
                             width=STD_POPULATION_V * 2,
                             height=STD_POPULATION_W * 2,
                             fill=False,
                             edgecolor=ELLIPSE_COLOR,
                             linestyle='dashed')
    uv_double_std_ellipse = Ellipse(xy=(AVERAGE_POPULATION_VELOCITY_U,
                                        AVERAGE_POPULATION_VELOCITY_V),
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_V * 4,
                                    fill=False,
                                    edgecolor=ELLIPSE_COLOR)
    uw_double_std_ellipse = Ellipse(xy=(AVERAGE_POPULATION_VELOCITY_U,
                                        AVERAGE_POPULATION_VELOCITY_W),
                                    width=STD_POPULATION_U * 4,
                                    height=STD_POPULATION_W * 4,
                                    fill=False,
                                    edgecolor=ELLIPSE_COLOR)
    vw_double_std_ellipse = Ellipse(xy=(AVERAGE_POPULATION_VELOCITY_V,
                                        AVERAGE_POPULATION_VELOCITY_W),
                                    width=STD_POPULATION_V * 4,
                                    height=STD_POPULATION_W * 4,
                                    fill=False,
                                    edgecolor=ELLIPSE_COLOR)

    subplot_uv.add_artist(uv_std_ellipse)
    subplot_uw.add_artist(uw_std_ellipse)
    subplot_vw.add_artist(vw_std_ellipse)
    subplot_uv.add_artist(uv_double_std_ellipse)
    subplot_uw.add_artist(uw_double_std_ellipse)
    subplot_vw.add_artist(vw_double_std_ellipse)

    # TODO: why does this apply minorticks only to the last subplot?
    plt.minorticks_on()

    subplot_uv.xaxis.set_ticks_position('both')
    subplot_uv.yaxis.set_ticks_position('both')
    subplot_uw.xaxis.set_ticks_position('both')
    subplot_uw.yaxis.set_ticks_position('both')
    subplot_vw.xaxis.set_ticks_position('both')
    subplot_vw.yaxis.set_ticks_position('both')

    subplot_uv.set_aspect(DESIRED_DIMENSIONS_RATIO
                          / subplot_uv.get_data_ratio())
    subplot_uw.set_aspect(DESIRED_DIMENSIONS_RATIO
                          / subplot_uw.get_data_ratio())
    subplot_vw.set_aspect(DESIRED_DIMENSIONS_RATIO
                          / subplot_vw.get_data_ratio())

    figure.subplots_adjust(hspace=SUBPLOTS_SPACING)

    plt.savefig(FILENAME)


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


def fetch_all_cloud_points(*,
                           session: Session):
    query = (Cloud.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [Cloud(**record)
            for record in records]


from typing import Tuple

from sqlalchemy.orm.session import Session
import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
from alcor.services.data_access import fetch_all_cloud_points
from alcor.services.data_access import (
    fetch_all_lepine_case_uv_cloud_points,
    fetch_all_lepine_case_uw_cloud_points,
    fetch_all_lepine_case_vw_cloud_points)

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.axes import Axes

FILENAME = 'velocity_clouds.ps'

# TODO: figure out how to work with sizes
FIGURE_SIZE = (8, 12)
DESIRED_DIMENSIONS_RATIO = 10 / 13

SUBPLOTS_SPACING = 0.25

# TODO: make dict for all labels and limits in all plotting modules
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


def plot(session: Session) -> None:
    # TODO: Implement other fetching functions
    cloud_points = fetch_all_cloud_points(session=session)

    velocities_u = [star.velocity_u
                    for star in cloud_points]
    velocities_v = [star.velocity_v
                    for star in cloud_points]
    velocities_w = [star.velocity_w
                    for star in cloud_points]

    figure, (subplot_uv,
             subplot_uw,
             subplot_vw) = plt.subplots(nrows=3,
                                        figsize=FIGURE_SIZE)

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

    plot_ellipses(subplot=subplot_uv,
                  xy=(AVERAGE_POPULATION_VELOCITY_U,
                      AVERAGE_POPULATION_VELOCITY_V),
                  x_std=STD_POPULATION_U,
                  y_std=STD_POPULATION_V)
    plot_ellipses(subplot=subplot_uw,
                  xy=(AVERAGE_POPULATION_VELOCITY_U,
                      AVERAGE_POPULATION_VELOCITY_W),
                  x_std=STD_POPULATION_U,
                  y_std=STD_POPULATION_W)
    plot_ellipses(subplot=subplot_vw,
                  xy=(AVERAGE_POPULATION_VELOCITY_V,
                      AVERAGE_POPULATION_VELOCITY_W),
                  x_std=STD_POPULATION_V,
                  y_std=STD_POPULATION_W)

    subplot_uv.minorticks_on()
    subplot_uw.minorticks_on()
    subplot_vw.minorticks_on()

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
    # TODO: Implement getting last points by time(ok?)
    uv_points = fetch_all_lepine_case_uv_cloud_points(session=session)
    uw_points = fetch_all_lepine_case_uw_cloud_points(session=session)
    vw_points = fetch_all_lepine_case_vw_cloud_points(session=session)

    uv_cloud_velocities_u = [star.velocity_u
                             for star in uv_points]
    uv_cloud_velocities_v = [star.velocity_v
                             for star in uv_points]
    uw_cloud_velocities_u = [star.velocity_u
                             for star in uw_points]
    uw_cloud_velocities_w = [star.velocity_w
                             for star in uw_points]
    vw_cloud_velocities_v = [star.velocity_v
                             for star in vw_points]
    vw_cloud_velocities_w = [star.velocity_w
                             for star in vw_points]

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

    subplot_uv.scatter(x=uv_cloud_velocities_u,
                       y=uv_cloud_velocities_v,
                       color=CLOUD_COLOR,
                       s=POINT_SIZE)
    subplot_uw.scatter(x=uw_cloud_velocities_u,
                       y=uw_cloud_velocities_w,
                       color=CLOUD_COLOR,
                       s=POINT_SIZE)
    subplot_vw.scatter(x=vw_cloud_velocities_v,
                       y=vw_cloud_velocities_w,
                       color=CLOUD_COLOR,
                       s=POINT_SIZE)

    plot_ellipses(subplot=subplot_uv,
                  xy=(AVERAGE_POPULATION_VELOCITY_U,
                      AVERAGE_POPULATION_VELOCITY_V),
                  x_std=STD_POPULATION_U,
                  y_std=STD_POPULATION_V)
    plot_ellipses(subplot=subplot_uw,
                  xy=(AVERAGE_POPULATION_VELOCITY_U,
                      AVERAGE_POPULATION_VELOCITY_W),
                  x_std=STD_POPULATION_U,
                  y_std=STD_POPULATION_W)
    plot_ellipses(subplot=subplot_vw,
                  xy=(AVERAGE_POPULATION_VELOCITY_V,
                      AVERAGE_POPULATION_VELOCITY_W),
                  x_std=STD_POPULATION_V,
                  y_std=STD_POPULATION_W)

    subplot_uv.minorticks_on()
    subplot_uw.minorticks_on()
    subplot_vw.minorticks_on()

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


def plot_ellipses(subplot: Axes,
                  xy: Tuple[float, float],
                  x_std: float,
                  y_std: float) -> None:
    std_ellipse = Ellipse(xy=xy,
                          width=x_std * 2,
                          height=y_std * 2,
                          fill=False,
                          edgecolor=ELLIPSE_COLOR,
                          linestyle='dashed')
    double_std_ellipse = Ellipse(xy=xy,
                                 width=x_std * 4,
                                 height=y_std * 4,
                                 fill=False,
                                 edgecolor=ELLIPSE_COLOR)

    subplot.add_artist(std_ellipse)
    subplot.add_artist(double_std_ellipse)

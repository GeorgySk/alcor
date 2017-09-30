from typing import (Tuple,
                    List)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.axes import Axes

from alcor.models import Star

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


def plot(stars: List[Star],
         filename: str = 'velocity_clouds.ps',
         figure_size: Tuple[float, float] = (8, 12),
         spacing: float = 0.25,
         u_label: str = '$U(km/s)$',
         v_label: str = '$V(km/s)$',
         w_label: str = '$W(km/s)$',
         u_limits: Tuple[float, float] = (-150, 150),
         v_limits: Tuple[float, float] = (-150, 150),
         w_limits: Tuple[float, float] = (-150, 150)) -> None:
    figure, (uv_subplot,
             uw_subplot,
             vw_subplot) = plt.subplots(nrows=3,
                                        figsize=figure_size)

    u_velocities = [star.u_velocity
                    for star in stars]
    v_velocities = [star.v_velocity
                    for star in stars]
    w_velocities = [star.w_velocity
                    for star in stars]

    draw_subplot(subplot=uv_subplot,
                 xlabel=u_label,
                 ylabel=v_label,
                 xlim=u_limits,
                 ylim=v_limits,
                 x=u_velocities,
                 y=v_velocities,
                 x_avg=AVERAGE_POPULATION_VELOCITY_U,
                 y_avg=AVERAGE_POPULATION_VELOCITY_V,
                 x_std=STD_POPULATION_U,
                 y_std=STD_POPULATION_V)
    draw_subplot(subplot=uw_subplot,
                 xlabel=u_label,
                 ylabel=w_label,
                 xlim=u_limits,
                 ylim=w_limits,
                 x=u_velocities,
                 y=w_velocities,
                 x_avg=AVERAGE_POPULATION_VELOCITY_U,
                 y_avg=AVERAGE_POPULATION_VELOCITY_W,
                 x_std=STD_POPULATION_U,
                 y_std=STD_POPULATION_W)
    draw_subplot(subplot=vw_subplot,
                 xlabel=v_label,
                 ylabel=w_label,
                 xlim=v_limits,
                 ylim=w_limits,
                 x=v_velocities,
                 y=w_velocities,
                 x_avg=AVERAGE_POPULATION_VELOCITY_V,
                 y_avg=AVERAGE_POPULATION_VELOCITY_W,
                 x_std=STD_POPULATION_V,
                 y_std=STD_POPULATION_W)

    figure.subplots_adjust(hspace=spacing)
    plt.savefig(filename)


def plot_lepine_case(stars: List[Star],
                     filename: str = 'velocity_clouds.ps',
                     figure_size: Tuple[float, float] = (8, 12),
                     spacing: float = 0.25,
                     u_label: str = '$U(km/s)$',
                     v_label: str = '$V(km/s)$',
                     w_label: str = '$W(km/s)$',
                     u_limits: Tuple[float, float] = (-150, 150),
                     v_limits: Tuple[float, float] = (-150, 150),
                     w_limits: Tuple[float, float] = (-150, 150)):
    figure, (uv_subplot,
             uw_subplot,
             vw_subplot) = plt.subplots(nrows=3,
                                        figsize=figure_size)

    uv_cloud_u_velocities = []
    uv_cloud_v_velocities = []
    uw_cloud_u_velocities = []
    uw_cloud_w_velocities = []
    vw_cloud_v_velocities = []
    vw_cloud_w_velocities = []

    for star in stars:
        max_coordinates_modulus = star.max_coordinates_modulus
        if abs(star.x_coordinate) == max_coordinates_modulus:
            vw_cloud_v_velocities.append(star.v_velocity)
            vw_cloud_w_velocities.append(star.w_velocity)
        elif abs(star.y_coordinate) == max_coordinates_modulus:
            uw_cloud_u_velocities.append(star.u_velocity)
            uw_cloud_w_velocities.append(star.w_velocity)
        else:
            uv_cloud_u_velocities.append(star.u_velocity)
            uv_cloud_v_velocities.append(star.v_velocity)

    draw_subplot(subplot=uv_subplot,
                 xlabel=u_label,
                 ylabel=v_label,
                 xlim=u_limits,
                 ylim=v_limits,
                 x=uv_cloud_u_velocities,
                 y=uv_cloud_v_velocities,
                 x_avg=AVERAGE_POPULATION_VELOCITY_U,
                 y_avg=AVERAGE_POPULATION_VELOCITY_V,
                 x_std=STD_POPULATION_U,
                 y_std=STD_POPULATION_V)
    draw_subplot(subplot=uw_subplot,
                 xlabel=u_label,
                 ylabel=w_label,
                 xlim=u_limits,
                 ylim=w_limits,
                 x=uw_cloud_u_velocities,
                 y=uw_cloud_w_velocities,
                 x_avg=AVERAGE_POPULATION_VELOCITY_U,
                 y_avg=AVERAGE_POPULATION_VELOCITY_W,
                 x_std=STD_POPULATION_U,
                 y_std=STD_POPULATION_W)
    draw_subplot(subplot=vw_subplot,
                 xlabel=v_label,
                 ylabel=w_label,
                 xlim=v_limits,
                 ylim=w_limits,
                 x=vw_cloud_v_velocities,
                 y=vw_cloud_w_velocities,
                 x_avg=AVERAGE_POPULATION_VELOCITY_V,
                 y_avg=AVERAGE_POPULATION_VELOCITY_W,
                 x_std=STD_POPULATION_V,
                 y_std=STD_POPULATION_W)

    figure.subplots_adjust(hspace=spacing)
    plt.savefig(filename)


def draw_subplot(*,
                 subplot: Axes,
                 xlabel: str,
                 ylabel: str,
                 xlim: Tuple[float, float],
                 ylim: Tuple[float, float],
                 x: List[float],
                 y: List[float],
                 cloud_color: str = 'k',
                 point_size: float = 0.5,
                 x_avg: float,
                 y_avg: float,
                 x_std: float,
                 y_std: float,
                 ratio: float = 10 / 13) -> None:
    subplot.set(xlabel=xlabel,
                ylabel=ylabel,
                xlim=xlim,
                ylim=ylim)
    subplot.scatter(x=x,
                    y=y,
                    color=cloud_color,
                    s=point_size)
    plot_ellipses(subplot=subplot,
                  x_avg=x_avg,
                  y_avg=y_avg,
                  x_std=x_std,
                  y_std=y_std)
    subplot.minorticks_on()
    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')
    subplot.set_aspect(ratio / subplot.get_data_ratio())


def plot_ellipses(subplot: Axes,
                  x_avg: float,
                  y_avg: float,
                  x_std: float,
                  y_std: float,
                  ellipse_color: str = 'b') -> None:
    std_ellipse = Ellipse(xy=(x_avg, y_avg),
                          width=x_std * 2,
                          height=y_std * 2,
                          fill=False,
                          edgecolor=ellipse_color,
                          linestyle='dashed')
    double_std_ellipse = Ellipse(xy=(x_avg, y_avg),
                                 width=x_std * 4,
                                 height=y_std * 4,
                                 fill=False,
                                 edgecolor=ellipse_color)

    subplot.add_artist(std_ellipse)
    subplot.add_artist(double_std_ellipse)

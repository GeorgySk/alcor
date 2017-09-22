import os
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
from sqlalchemy.orm.session import Session
import pandas as pd

from alcor.services.data_access import fetch_all
from alcor.models.velocities.clouds import (Cloud,
                                            LepineCaseUVCloud,
                                            LepineCaseUWCloud,
                                            LepineCaseVWCloud)

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


def plot(session: Session,
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

    # TODO: Implement other fetching functions
    cloud_points = fetch_all(Cloud,
                             session=session)

    u_velocities = [star.u_velocity
                    for star in cloud_points]
    v_velocities = [star.v_velocity
                    for star in cloud_points]
    w_velocities = [star.w_velocity
                    for star in cloud_points]

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


def plot_lepine_case(session: Session,
                     filename: str = 'velocity_clouds.ps',
                     figure_size: Tuple[float, float] = (8, 12),
                     spacing: float = 0.25,
                     u_label: str = '$U(km/s)$',
                     v_label: str = '$V(km/s)$',
                     w_label: str = '$W(km/s)$',
                     u_limits: Tuple[float, float] = (-150, 150),
                     v_limits: Tuple[float, float] = (-150, 150),
                     w_limits: Tuple[float, float] = (-150, 150),
                     uv_observational_path: str = 'observational_data/uv.dat',
                     uw_observational_path: str = 'observational_data/uw.dat',
                     vw_observational_path: str = 'observational_data/vw.dat'
                     ) -> None:
    base_dir = os.path.dirname(__file__)
    observational_clouds_paths = dict(
            uv=os.path.join(base_dir, uv_observational_path),
            uw=os.path.join(base_dir, uw_observational_path),
            vw=os.path.join(base_dir, vw_observational_path))

    uv_observational_cloud_df = pd.read_csv(observational_clouds_paths['uv'],
                                            delimiter=' ',
                                            skipinitialspace=True,
                                            header=None)
    uw_observational_cloud_df = pd.read_csv(observational_clouds_paths['uw'],
                                            delimiter=' ',
                                            skipinitialspace=True,
                                            header=None)
    vw_observational_cloud_df = pd.read_csv(observational_clouds_paths['vw'],
                                            delimiter=' ',
                                            skipinitialspace=True,
                                            header=None)

    # TODO: Add other fetching options
    uv_points = fetch_all(LepineCaseUVCloud,
                          session=session)
    uw_points = fetch_all(LepineCaseUWCloud,
                          session=session)
    vw_points = fetch_all(LepineCaseVWCloud,
                          session=session)

    uv_cloud_u_velocities = [star.u_velocity
                             for star in uv_points]
    uv_cloud_v_velocities = [star.v_velocity
                             for star in uv_points]

    uw_cloud_u_velocities = [star.u_velocity
                             for star in uw_points]
    uw_cloud_w_velocities = [star.w_velocity
                             for star in uw_points]

    vw_cloud_v_velocities = [star.v_velocity
                             for star in vw_points]
    vw_cloud_w_velocities = [star.w_velocity
                             for star in vw_points]

    figure, (uv_subplot,
             uw_subplot,
             vw_subplot) = plt.subplots(nrows=3,
                                        figsize=figure_size)

    draw_subplot(subplot=uv_subplot,
                 xlabel=u_label,
                 ylabel=v_label,
                 xlim=u_limits,
                 ylim=v_limits,
                 x=uv_cloud_u_velocities,
                 y=uv_cloud_v_velocities,
                 x_obs=uv_observational_cloud_df[0],
                 y_obs=uv_observational_cloud_df[1],
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
                 x_obs=uw_observational_cloud_df[0],
                 y_obs=uw_observational_cloud_df[1],
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
                 x_obs=vw_observational_cloud_df[0],
                 y_obs=vw_observational_cloud_df[1],
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
                 x_obs: pd.core.series.Series = None,
                 y_obs: pd.core.series.Series = None,
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
    subplot.scatter(x=x_obs,
                    y=y_obs,
                    color='r',
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

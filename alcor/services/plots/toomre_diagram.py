import enum
import logging
from math import sqrt
from typing import (Tuple,
                    List)

import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from sqlalchemy.orm.session import Session

from alcor.models import (Star,
                          GalacticDiskType)
from alcor.services.common import PECULIAR_SOLAR_VELOCITY_V
from alcor.services.data_access.service import fetch_random

logger = logging.getLogger(__name__)


def plot(session: Session,
         filename: str = 'toomre_diagram.ps',
         figure_size: Tuple[float, float] = (8, 8),
         ratio: float = 10 / 13,
         xlabel: str = '$V(km/s)$',
         ylabel: str = '$\sqrt{U^2+W^2}(km/s)$',
         thin_disk_color: str = 'r',
         thick_disk_color: str = 'b',
         desired_stars_count: int = 10_000) -> None:
    figure, subplot = plt.subplots(figsize=figure_size)

    # TODO: add other fetching options
    stars = fetch_random(Star,
                         limit=desired_stars_count,
                         session=session)

    # TODO: add choosing frame: relative to Sun/LSR. Now it's rel. to LSR
    plot_stars_by_galactic_disk_type(subplot=subplot,
                                     stars=stars,
                                     galactic_disk_type=GalacticDiskType.thin,
                                     color=thin_disk_color)
    plot_stars_by_galactic_disk_type(subplot=subplot,
                                     stars=stars,
                                     galactic_disk_type=GalacticDiskType.thick,
                                     color=thick_disk_color)

    # TODO: add sliders
    subplot.set(xlabel=xlabel,
                ylabel=ylabel)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    plt.savefig(filename)


def plot_stars_by_galactic_disk_type(*,
                                     subplot: Axes,
                                     stars: List[Star],
                                     galactic_disk_type: enum.Enum,
                                     color: str,
                                     point_size: float = 0.5) -> None:
    # TODO: how to work with Decimal type? If I leave it I get:
    # TypeError: Cannot cast array data from dtype('O') to dtype('float64')
    # according to the rule 'safe'
    v_velocities = [float(star.v_velocity)
                    + PECULIAR_SOLAR_VELOCITY_V
                    for star in stars
                    if star.galactic_disk_type == galactic_disk_type]
    uw_velocities_square_sums_square_root = [
        sqrt(float(star.u_velocity) ** 2 + float(star.w_velocity) ** 2)
        for star in stars
        if star.galactic_disk_type == galactic_disk_type]

    subplot.scatter(x=v_velocities,
                    y=uw_velocities_square_sums_square_root,
                    color=color,
                    s=point_size)

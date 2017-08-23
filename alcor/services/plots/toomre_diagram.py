import enum
import logging
from math import sqrt
from random import random
from typing import (Tuple,
                    List)

from matplotlib.axes import Axes
from sqlalchemy.orm.session import Session
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470

from alcor.models.star import (Star,
                               GalacticDiskEnum)
from alcor.services.data_access import fetch_all_stars
from alcor.services.restrictions import PECULIAR_SOLAR_VELOCITY_V

logger = logging.getLogger(__name__)

DESIRED_FINAL_SAMPLE_STARS_COUNT = 10_000


def plot(session: Session,
         filename: str = 'toomre_diagram.ps',
         figure_size: Tuple[float, float] = (8, 8),
         ratio: float = 10 / 13,
         xlabel: str = '$V(km/s)$',
         ylabel: str = '$\sqrt{U^2+W^2}(km/s)$',
         thin_disk_color: str = 'r',
         thick_disk_color: str = 'b') -> None:
    figure, subplot = plt.subplots(figsize=figure_size)

    # TODO: add other fetching options
    stars = fetch_all_stars(session=session)

    # TODO: or is it better to place this in fetch_.. function?
    choosing_probability = DESIRED_FINAL_SAMPLE_STARS_COUNT / len(stars)
    random_stars_sample = [star for star in stars
                           if random() < choosing_probability]

    # TODO: add choosing frame: relative to Sun/LSR. Now it's rel. to LSR
    plot_stars_by_disk(subplot=subplot,
                       stars=random_stars_sample,
                       galactic_disk=GalacticDiskEnum.thin,
                       color=thin_disk_color)
    plot_stars_by_disk(subplot=subplot,
                       stars=random_stars_sample,
                       galactic_disk=GalacticDiskEnum.thick,
                       color=thick_disk_color)

    # TODO: add sliders
    subplot.set(xlabel=xlabel,
                ylabel=ylabel)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    plt.savefig(filename)


def plot_stars_by_disk(*,
                       subplot: Axes,
                       stars: List[Star],
                       galactic_disk: enum.Enum,
                       color: str,
                       point_size: float = 0.5) -> None:
    # TODO: how to work with Decimal type? If I leave it I get:
    # TypeError: Cannot cast array data from dtype('O') to dtype('float64')
    # according to the rule 'safe'
    velocities_v = [float(star.velocity_v)
                    + PECULIAR_SOLAR_VELOCITY_V
                    for star in stars
                    if star.disk_belonging == galactic_disk]
    uw_velocities_square_sums_square_root = [
        sqrt(float(star.velocity_u) ** 2 + float(star.velocity_w) ** 2)
        for star in stars
        if star.disk_belonging == galactic_disk]

    subplot.scatter(x=velocities_v,
                    y=uw_velocities_square_sums_square_root,
                    color=color,
                    s=point_size)

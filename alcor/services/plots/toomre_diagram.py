import logging
from math import sqrt
from random import random
from typing import List

from matplotlib.axes import Axes
from sqlalchemy.orm.session import Session
import matplotlib

# More info at
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
from alcor.models.star import Star
from alcor.services.data_access import fetch_all_stars

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.services.restrictions import PECULIAR_SOLAR_VELOCITY_V

logger = logging.getLogger(__name__)

THIN_DISK_INDEX = 1
THICK_DISK_INDEX = 2

DESIRED_FINAL_SAMPLE_STARS_COUNT = 10_000

FILENAME = 'toomre_diagram.ps'

FIGURE_SIZE = (8, 8)
DESIRED_DIMENSIONS_RATIO = 10 / 13

X_LABEL = '$V(km/s)$'
Y_LABEL = '$\sqrt{U^2+W^2}(km/s)$'

THIN_DISK_CLOUD_COLOR = 'r'
THICK_DISK_CLOUD_COLOR = 'b'


def plot(session: Session) -> None:
    figure, subplot = plt.subplots(figsize=FIGURE_SIZE)

    # TODO: add other fetching options
    stars = fetch_all_stars(session=session)

    # TODO: or is it better to place this in fetch_.. function?
    choosing_probability = DESIRED_FINAL_SAMPLE_STARS_COUNT / len(stars)
    random_stars_sample = [star for star in stars
                           if random() < choosing_probability]

    # TODO: add choosing frame: relative to Sun/LSR. Now it's rel. to LSR
    plot_stars_by_disk(subplot=subplot,
                       stars=random_stars_sample,
                       disk_index=THIN_DISK_INDEX,
                       color=THIN_DISK_CLOUD_COLOR)
    plot_stars_by_disk(subplot=subplot,
                       stars=random_stars_sample,
                       disk_index=THICK_DISK_INDEX,
                       color=THICK_DISK_CLOUD_COLOR)

    # TODO: add sliders
    subplot.set(xlabel=X_LABEL,
                ylabel=Y_LABEL)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(DESIRED_DIMENSIONS_RATIO
                       / subplot.get_data_ratio())

    plt.savefig(FILENAME)


def plot_stars_by_disk(*,
                       subplot: Axes,
                       stars: List[Star],
                       disk_index: int,
                       color: str,
                       point_size: float = 0.5) -> None:
    # TODO: how to work with Decimal type? If I leave it I get:
    # TypeError: Cannot cast array data from dtype('O') to dtype('float64')
    # according to the rule 'safe'
    velocities_v = [float(star.velocity_v)
                    + PECULIAR_SOLAR_VELOCITY_V
                    for star in stars
                    if star.disk_belonging == disk_index]
    uw_velocities_square_sums_square_root = [
        sqrt(float(star.velocity_u) ** 2 + float(star.velocity_w) ** 2)
        for star in stars
        if star.disk_belonging == disk_index]

    subplot.scatter(x=velocities_v,
                    y=uw_velocities_square_sums_square_root,
                    color=color,
                    s=point_size)

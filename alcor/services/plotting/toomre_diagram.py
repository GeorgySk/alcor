import logging
from math import sqrt
from random import random

from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.models.star import Star
from alcor.services.data_access.reading import fetch


logger = logging.getLogger(__name__)


DESIRED_FINAL_SAMPLE_STARS_COUNT = 10_000

FILENAME = 'toomre_diagram.ps'

FIGURE_SIZE = (8, 8)
DESIRED_DIMENSIONS_RATIO = 10 / 13

X_LABEL = '$V(km/s)$'
Y_LABEL = '$\sqrt{U^2+W^2}(km/s)$'

# TODO: in what reference frame?
PECULIAR_SOLAR_VELOCITY_U = -11
PECULIAR_SOLAR_VELOCITY_V = 12
PECULIAR_SOLAR_VELOCITY_W = 7

THIN_DISK_CLOUD_COLOR = 'r'
THICK_DISK_CLOUD_COLOR = 'b'
POINT_SIZE = 0.5


def plot(*,
         session: Session) -> None:
    # TODO: Figure out what stars I should fetch (all/last group by time/last N
    # groups by time/selected by ID/marked by some flag(series of simulations))
    stars = fetch_all_stars(session=session)

    # TODO: or is it better to place this in fetch_.. function?
    choosing_probability = DESIRED_FINAL_SAMPLE_STARS_COUNT / len(stars)
    random_stars_sample = [star for star in stars
                           if random() < choosing_probability]

    # TODO: add choosing frame: relative to Sun/LSR. Now it's rel. to LSR
    # TODO: how to work with Decimal type? If I leave it I get:
    # TypeError: Cannot cast array data from dtype('O') to dtype('float64')
    # according to the rule 'safe'
    thin_disk_velocities_v = [float(star.velocity_v)
                              + PECULIAR_SOLAR_VELOCITY_V
                              for star in random_stars_sample
                              if star.disk_belonging == 1]
    thin_disk_uw_velocities_square_sums_square_root = [
        sqrt(float(star.velocity_u) ** 2 + float(star.velocity_w) ** 2)
        for star in random_stars_sample
        if star.disk_belonging == 1]

    thick_disk_velocities_v = [float(star.velocity_v)
                               + PECULIAR_SOLAR_VELOCITY_V
                               for star in random_stars_sample
                               if star.disk_belonging == 2]
    thick_disk_uw_velocities_square_sums_square_root = [
        sqrt(float(star.velocity_u) ** 2 + float(star.velocity_w) ** 2)
        for star in random_stars_sample
        if star.disk_belonging == 2]

    figure, subplot = plt.subplots(figsize=FIGURE_SIZE)

    # TODO: add sliders
    subplot.set(xlabel=X_LABEL,
                ylabel=Y_LABEL)

    subplot.scatter(x=thick_disk_velocities_v,
                    y=thick_disk_uw_velocities_square_sums_square_root,
                    color=THICK_DISK_CLOUD_COLOR,
                    s=POINT_SIZE)
    subplot.scatter(x=thin_disk_velocities_v,
                    y=thin_disk_uw_velocities_square_sums_square_root,
                    color=THIN_DISK_CLOUD_COLOR,
                    s=POINT_SIZE)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(DESIRED_DIMENSIONS_RATIO
                       / subplot.get_data_ratio())

    plt.savefig(FILENAME)


def fetch_all_stars(*,
                    session: Session):
    query = (Star.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [Star(**record)
            for record in records]

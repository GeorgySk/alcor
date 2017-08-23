import logging
from statistics import (mean,
                        stdev)
from typing import (Iterator,
                    List)

from alcor.models import (Group,
                          Star)
from alcor.models.velocities_vs_magnitudes import Bin
from alcor.types import StarsBinsType
from .utils import (BIN_SIZE,
                    BINS_COUNT,
                    DEFAULT_VELOCITY_STD,
                    MIN_BOLOMETRIC_MAGNITUDE,
                    get_stars_bin_index)

logger = logging.getLogger(__name__)


def bins(*,
         stars: List[Star],
         group: Group) -> Iterator[Bin]:
    group_id = group.id
    for index, stars_bin in enumerate(stars_bins(stars)):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (index + 0.5))
        avg_u_velocity = mean(star.u_velocity
                              for star in stars_bin)
        avg_v_velocity = mean(star.v_velocity
                              for star in stars_bin)
        avg_w_velocity = mean(star.w_velocity
                              for star in stars_bin)
        if len(stars_bin) > 1:
            u_velocity_std = stdev(star.u_velocity
                                   for star in stars_bin)
            v_velocity_std = stdev(star.v_velocity
                                   for star in stars_bin)
            w_velocity_std = stdev(star.w_velocity
                                   for star in stars_bin)
        else:
            u_velocity_std = DEFAULT_VELOCITY_STD
            v_velocity_std = DEFAULT_VELOCITY_STD
            w_velocity_std = DEFAULT_VELOCITY_STD
        yield Bin(group_id=group_id,
                  avg_magnitude=avg_magnitude,
                  avg_u_velocity=avg_u_velocity,
                  avg_v_velocity=avg_v_velocity,
                  avg_w_velocity=avg_w_velocity,
                  u_velocity_std=u_velocity_std,
                  v_velocity_std=v_velocity_std,
                  w_velocity_std=w_velocity_std)


def stars_bins(stars: List[Star]) -> StarsBinsType:
    res = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        index = get_stars_bin_index(star)
        if BINS_COUNT > index >= 0:
            res[index].append(star)
        else:
            logger.warning('Magnitude is out of bounds: '
                           f'{star.bolometric_magnitude}')
    return res

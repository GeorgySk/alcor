from statistics import (mean,
                        stdev)
from typing import (Iterator,
                    List)

from alcor.models import Star
from alcor.models.velocities_vs_magnitudes import Bin
from .utils import (STARS_BIN_SIZE,
                    DEFAULT_VELOCITY_STD,
                    MIN_BOLOMETRIC_MAGNITUDE,
                    pack_stars)


def bins(stars: List[Star]) -> Iterator[Bin]:
    stars_bins = pack_stars(stars)
    for index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + STARS_BIN_SIZE * (index + 0.5))
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
        yield Bin(avg_magnitude=avg_magnitude,
                  avg_u_velocity=avg_u_velocity,
                  avg_v_velocity=avg_v_velocity,
                  avg_w_velocity=avg_w_velocity,
                  u_velocity_std=u_velocity_std,
                  v_velocity_std=v_velocity_std,
                  w_velocity_std=w_velocity_std)

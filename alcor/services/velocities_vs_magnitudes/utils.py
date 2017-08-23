import logging
from statistics import (mean,
                        stdev)
from typing import (Union,
                    Iterable,
                    Iterator,
                    Tuple,
                    List)

from alcor.models import (Group,
                          Star)
from alcor.models.velocities_vs_magnitudes import (Bin,
                                                   LepineCaseUBin,
                                                   LepineCaseVBin,
                                                   LepineCaseWBin,
                                                   LepineCaseUCloud,
                                                   LepineCaseVCloud,
                                                   LepineCaseWCloud)
from alcor.types import (StarsBinsType,
                         RowType)

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 30.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)
DEFAULT_VELOCITY_STD = 100.

logger = logging.getLogger(__name__)


def stars_clouds(stars: List[Star],
                 group: Group) -> Iterator[Union[LepineCaseUCloud,
                                                 LepineCaseVCloud,
                                                 LepineCaseWCloud]]:
    group_id = group.id
    for star in stars:
        max_coordinates_modulus = star.max_coordinates_modulus
        bolometric_magnitude = star.bolometric_magnitude

        if abs(star.x_coordinate) == max_coordinates_modulus:
            yield LepineCaseVCloud(group_id=group_id,
                                   v_velocity=star.v_velocity,
                                   bolometric_magnitude=bolometric_magnitude)
            yield LepineCaseWCloud(group_id=group_id,
                                   w_velocity=star.w_velocity,
                                   bolometric_magnitude=bolometric_magnitude)
        elif abs(star.y_coordinate) == max_coordinates_modulus:
            yield LepineCaseUCloud(group_id=group_id,
                                   u_velocity=star.u_velocity,
                                   bolometric_magnitude=bolometric_magnitude)
            yield LepineCaseWCloud(group_id=group_id,
                                   w_velocity=star.w_velocity,
                                   bolometric_magnitude=bolometric_magnitude)
        else:
            yield LepineCaseUCloud(group_id=group_id,
                                   u_velocity=star.u_velocity,
                                   bolometric_magnitude=bolometric_magnitude)
            yield LepineCaseVCloud(group_id=group_id,
                                   v_velocity=star.v_velocity,
                                   bolometric_magnitude=bolometric_magnitude)


def generate_u_bins(*,
                    stars_bins: StarsBinsType,
                    group: Group) -> Iterable[LepineCaseUBin]:
    for index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (index + 0.5))
        avg_u_velocity = mean(star.u_velocity
                              for star in stars_bin)
        if len(stars_bin) == 1:
            u_velocity_std = DEFAULT_VELOCITY_STD
        else:
            u_velocity_std = stdev(star.u_velocity
                                   for star in stars_bin)
        yield LepineCaseUBin(group_id=group.id,
                             avg_magnitude=avg_magnitude,
                             avg_u_velocity=avg_u_velocity,
                             u_velocity_std=u_velocity_std)


def generate_v_bins(*,
                    stars_bins: StarsBinsType,
                    group: Group) -> Iterable[LepineCaseVBin]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index + 0.5))
        avg_v_velocity = mean(star.v_velocity
                              for star in stars_bin)
        if len(stars_bin) == 1:
            v_velocity_std = DEFAULT_VELOCITY_STD
        else:
            v_velocity_std = stdev(star.v_velocity
                                   for star in stars_bin)
        yield LepineCaseVBin(group_id=group.id,
                             avg_magnitude=avg_magnitude,
                             avg_v_velocity=avg_v_velocity,
                             v_velocity_std=v_velocity_std)


def generate_w_bins(*,
                    stars_bins: StarsBinsType,
                    group: Group) -> Iterable[RowType]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index + 0.5))
        avg_w_velocity = mean(star.w_velocity
                              for star in stars_bin)
        if len(stars_bin) == 1:
            w_velocity_std = DEFAULT_VELOCITY_STD
        else:
            w_velocity_std = stdev(star.w_velocity
                                   for star in stars_bin)
        yield LepineCaseWBin(group_id=group.id,
                             avg_magnitude=avg_magnitude,
                             avg_w_velocity=avg_w_velocity,
                             w_velocity_std=w_velocity_std)


def lepine_stars_bins(stars: List[Star]) -> Tuple[StarsBinsType,
                                                  StarsBinsType,
                                                  StarsBinsType]:
    u_stars_bins = [[] for _ in range(BINS_COUNT)]
    v_stars_bins = [[] for _ in range(BINS_COUNT)]
    w_stars_bins = [[] for _ in range(BINS_COUNT)]

    for star in stars:
        index = get_stars_bin_index(star)

        max_coordinates_modulus = star.max_coordinates_modulus

        if abs(star.x_coordinate) == max_coordinates_modulus:
            v_stars_bins[index].append(star)
            w_stars_bins[index].append(star)
        elif abs(star.y_coordinate) == max_coordinates_modulus:
            u_stars_bins[index].append(star)
            w_stars_bins[index].append(star)
        else:
            u_stars_bins[index].append(star)
            v_stars_bins[index].append(star)

    return u_stars_bins, v_stars_bins, w_stars_bins


def get_stars_bin_index(star: Star) -> int:
    return int((float(star.bolometric_magnitude)
                - MIN_BOLOMETRIC_MAGNITUDE) / BIN_SIZE)


def raw_stars_bins(stars: List[Star]) -> StarsBinsType:
    res = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        index = get_stars_bin_index(star)
        if BINS_COUNT > index >= 0:
            res[index].append(star)
        else:
            logger.warning(f'Magnitude is out of bounds: '
                           f'{star.bolometric_magnitude}')
    return res


def generate_bins(*,
                  stars_bins: StarsBinsType,
                  group: Group) -> Iterable[Bin]:
    for index, stars_bin in enumerate(stars_bins):
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
        yield Bin(group_id=group.id,
                  avg_magnitude=avg_magnitude,
                  avg_u_velocity=avg_u_velocity,
                  avg_v_velocity=avg_v_velocity,
                  avg_w_velocity=avg_w_velocity,
                  u_velocity_std=u_velocity_std,
                  v_velocity_std=v_velocity_std,
                  w_velocity_std=w_velocity_std)

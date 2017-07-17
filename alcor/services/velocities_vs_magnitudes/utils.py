from math import ceil
from statistics import (mean,
                        stdev)
from typing import (Iterable,
                    Tuple,
                    List)

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.cassandra_models.velocities_vs_magnitudes import (CBin,
                                                             CLepineCaseUBin,
                                                             CLepineCaseVBin,
                                                             CLepineCaseUCloud,
                                                             CLepineCaseVCloud,
                                                             CLepineCaseWCloud)
from alcor.types import (RowType,
                         StarsBinsType)

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)
DEFAULT_VELOCITY_STD = 100.


def generate_clouds(stars: List[CStar],
                    group: CGroup,
                    u_cls,
                    v_cls,
                    w_cls) -> Tuple[List[CLepineCaseUCloud],
                                    List[CLepineCaseVCloud],
                                    List[CLepineCaseWCloud]]:
    u_clouds = []
    v_clouds = []
    w_clouds = []

    for star in stars:
        highest_coordinate = max(star.coordinate_x,
                                 star.coordinate_y,
                                 star.coordinate_z)
        if star.coordinate_x == highest_coordinate:
            u_clouds.append(u_cls(
                group_id=group.id,
                bolometric_magnitude=star.bolometric_magnitude,
                velocity_u=star.velocity_u))
        elif star.coordinate_y == highest_coordinate:
            v_clouds.append(v_cls(
                group_id=group.id,
                bolometric_magnitude=star.bolometric_magnitude,
                velocity_v=star.velocity_v))
        else:
            w_clouds.append(w_cls(
                group_id=group.id,
                bolometric_magnitude=star.bolometric_magnitude,
                velocity_w=star.velocity_w))

    return u_clouds, v_clouds, w_clouds


def generate_u_bins(*,
                    stars_bins: StarsBinsType,
                    group: CGroup,
                    cls) -> Iterable[CLepineCaseUBin]:
    non_empty_stars_bins = filter(None, stars_bins)
    for stars_bin_index, stars_bin in enumerate(non_empty_stars_bins):
        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index - 0.5))
        avg_velocity_u = mean(star.velocity_u
                              for star in stars_bin)
        if len(stars_bin) == 1:
            velocity_u_std = DEFAULT_VELOCITY_STD
        else:
            velocity_u_std = stdev(star.velocity_u
                                   for star in stars_bin)
        yield cls(group_id=group.id,
                  avg_magnitude=avg_magnitude,
                  avg_velocity_u=avg_velocity_u,
                  velocity_u_std=velocity_u_std)


def generate_v_bins(*,
                    stars_bins: StarsBinsType,
                    group: CGroup,
                    cls) -> Iterable[CLepineCaseVBin]:
    non_empty_stars_bins = filter(None, stars_bins)
    for stars_bin_index, stars_bin in enumerate(non_empty_stars_bins):
        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index - 0.5))
        avg_velocity_v = mean(star.velocity_v
                              for star in stars_bin)
        if len(stars_bin) == 1:
            velocity_v_std = DEFAULT_VELOCITY_STD
        else:
            velocity_v_std = stdev(star.velocity_v
                                   for star in stars_bin)
        yield cls(group_id=group.id,
                  avg_magnitude=avg_magnitude,
                  avg_velocity_v=avg_velocity_v,
                  velocity_v_std=velocity_v_std)


def generate_w_bins(*,
                    stars_bins: StarsBinsType,
                    group: CGroup,
                    cls) -> Iterable[RowType]:
    non_empty_stars_bins = filter(None, stars_bins)
    for stars_bin_index, stars_bin in enumerate(non_empty_stars_bins):
        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index - 0.5))
        avg_velocity_w = mean(star.velocity_w
                              for star in stars_bin)
        if len(stars_bin) == 1:
            velocity_w_std = DEFAULT_VELOCITY_STD
        else:
            velocity_w_std = stdev(star.velocity_w
                                   for star in stars_bin)
        yield cls(group_id=group.id,
                  avg_magnitude=avg_magnitude,
                  avg_velocity_w=avg_velocity_w,
                  velocity_w_std=velocity_w_std)


def generate_stars_bins(stars: List[CStar]) -> Tuple[StarsBinsType,
                                                     StarsBinsType,
                                                     StarsBinsType]:
    u_stars_bins = [[] for _ in range(BINS_COUNT)]
    v_stars_bins = [[] for _ in range(BINS_COUNT)]
    w_stars_bins = [[] for _ in range(BINS_COUNT)]

    for star in stars:
        index = get_stars_bin_index(star)

        highest_coordinate = max(star.coordinate_x,
                                 star.coordinate_y,
                                 star.coordinate_z)
        if star.coordinate_x == highest_coordinate:
            v_stars_bins[index].append(star)
            w_stars_bins[index].append(star)
        elif star.coordinate_y == highest_coordinate:
            u_stars_bins[index].append(star)
            w_stars_bins[index].append(star)
        else:
            u_stars_bins[index].append(star)
            v_stars_bins[index].append(star)

    return u_stars_bins, v_stars_bins, w_stars_bins


def get_stars_bin_index(star: CStar) -> int:
    return int(ceil((float(star.bolometric_magnitude)
                     - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))


def raw_stars_bins(stars: List[CStar]) -> StarsBinsType:
    res = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        index = get_stars_bin_index(star)
        res[index].append(star)
    return res


def generate_bins(*,
                  stars_bins: StarsBinsType,
                  group: CGroup,
                  cls) -> Iterable[CBin]:
    non_empty_stars_bins = filter(None, stars_bins)
    for index, stars_bin in enumerate(non_empty_stars_bins):
        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (index - 0.5))
        avg_velocity_u = mean(star.velocity_u
                              for star in stars_bin)
        avg_velocity_v = mean(star.velocity_v
                              for star in stars_bin)
        avg_velocity_w = mean(star.velocity_w
                              for star in stars_bin)
        if len(stars_bin) > 1:
            velocity_u_std = stdev(star.velocity_u
                                   for star in stars_bin)
            velocity_v_std = stdev(star.velocity_v
                                   for star in stars_bin)
            velocity_w_std = stdev(star.velocity_w
                                   for star in stars_bin)
        else:
            velocity_u_std = DEFAULT_VELOCITY_STD
            velocity_v_std = DEFAULT_VELOCITY_STD
            velocity_w_std = DEFAULT_VELOCITY_STD
        yield cls(group_id=group.id,
                  avg_magnitude=avg_magnitude,
                  avg_velocity_u=avg_velocity_u,
                  avg_velocity_v=avg_velocity_v,
                  avg_velocity_w=avg_velocity_w,
                  velocity_u_std=velocity_u_std,
                  velocity_v_std=velocity_v_std,
                  velocity_w_std=velocity_w_std)

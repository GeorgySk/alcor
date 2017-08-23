import logging
from statistics import (mean,
                        stdev)
from typing import (Union,
                    Iterable,
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


def generate_clouds(stars: List[Star],
                    group: Group) -> List[Union[LepineCaseUCloud,
                                                LepineCaseVCloud,
                                                LepineCaseWCloud]]:
    u_clouds = []
    v_clouds = []
    w_clouds = []
    for star in stars:
        max_coordinates_modulus = max(abs(star.coordinate_x),
                                      abs(star.coordinate_y),
                                      abs(star.coordinate_z))
        if abs(star.coordinate_x) == max_coordinates_modulus:
            v_clouds.append(
                LepineCaseVCloud(group_id=group.id,
                                 velocity_v=star.velocity_v,
                                 bolometric_magnitude=star.bolometric_magnitude))
            w_clouds.append(
                LepineCaseWCloud(group_id=group.id,
                                 velocity_w=star.velocity_w,
                                 bolometric_magnitude=star.bolometric_magnitude))
        elif abs(star.coordinate_y) == max_coordinates_modulus:
            u_clouds.append(
                LepineCaseUCloud(group_id=group.id,
                                 velocity_u=star.velocity_u,
                                 bolometric_magnitude=star.bolometric_magnitude))
            w_clouds.append(
                LepineCaseWCloud(group_id=group.id,
                                 velocity_w=star.velocity_w,
                                 bolometric_magnitude=star.bolometric_magnitude))
        else:
            u_clouds.append(
                LepineCaseUCloud(group_id=group.id,
                                 velocity_u=star.velocity_u,
                                 bolometric_magnitude=star.bolometric_magnitude))
            v_clouds.append(
                LepineCaseVCloud(group_id=group.id,
                                 velocity_v=star.velocity_v,
                                 bolometric_magnitude=star.bolometric_magnitude))
    clouds = []
    clouds.extend(u_clouds)
    clouds.extend(v_clouds)
    clouds.extend(w_clouds)

    return clouds


def generate_u_bins(*,
                    stars_bins: StarsBinsType,
                    group: Group) -> Iterable[LepineCaseUBin]:
    for index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (index + 0.5))
        avg_velocity_u = mean(star.velocity_u
                              for star in stars_bin)
        if len(stars_bin) == 1:
            velocity_u_std = DEFAULT_VELOCITY_STD
        else:
            velocity_u_std = stdev(star.velocity_u
                                   for star in stars_bin)
        yield LepineCaseUBin(group_id=group.id,
                             avg_magnitude=avg_magnitude,
                             avg_velocity_u=avg_velocity_u,
                             velocity_u_std=velocity_u_std)


def generate_v_bins(*,
                    stars_bins: StarsBinsType,
                    group: Group) -> Iterable[LepineCaseVBin]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index + 0.5))
        avg_velocity_v = mean(star.velocity_v
                              for star in stars_bin)
        if len(stars_bin) == 1:
            velocity_v_std = DEFAULT_VELOCITY_STD
        else:
            velocity_v_std = stdev(star.velocity_v
                                   for star in stars_bin)
        yield LepineCaseVBin(group_id=group.id,
                             avg_magnitude=avg_magnitude,
                             avg_velocity_v=avg_velocity_v,
                             velocity_v_std=velocity_v_std)


def generate_w_bins(*,
                    stars_bins: StarsBinsType,
                    group: Group) -> Iterable[RowType]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + BIN_SIZE * (stars_bin_index + 0.5))
        avg_velocity_w = mean(star.velocity_w
                              for star in stars_bin)
        if len(stars_bin) == 1:
            velocity_w_std = DEFAULT_VELOCITY_STD
        else:
            velocity_w_std = stdev(star.velocity_w
                                   for star in stars_bin)
        yield LepineCaseWBin(group_id=group.id,
                             avg_magnitude=avg_magnitude,
                             avg_velocity_w=avg_velocity_w,
                             velocity_w_std=velocity_w_std)


def lepine_stars_bins(stars: List[Star]) -> Tuple[StarsBinsType,
                                                  StarsBinsType,
                                                  StarsBinsType]:
    u_stars_bins = [[] for _ in range(BINS_COUNT)]
    v_stars_bins = [[] for _ in range(BINS_COUNT)]
    w_stars_bins = [[] for _ in range(BINS_COUNT)]

    for star in stars:
        index = get_stars_bin_index(star)

        max_coordinates_modulus = max(abs(star.coordinate_x),
                                      abs(star.coordinate_y),
                                      abs(star.coordinate_z))
        if abs(star.coordinate_x) == max_coordinates_modulus:
            v_stars_bins[index].append(star)
            w_stars_bins[index].append(star)
        elif abs(star.coordinate_y) == max_coordinates_modulus:
            u_stars_bins[index].append(star)
            w_stars_bins[index].append(star)
        else:
            u_stars_bins[index].append(star)
            v_stars_bins[index].append(star)

    return u_stars_bins, v_stars_bins, w_stars_bins


def get_stars_bin_index(star: Star) -> int:
    return int((float(star.bolometric_magnitude)
                - MIN_BOLOMETRIC_MAGNITUDE)/BIN_SIZE)


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
        yield Bin(group_id=group.id,
                  avg_magnitude=avg_magnitude,
                  avg_velocity_u=avg_velocity_u,
                  avg_velocity_v=avg_velocity_v,
                  avg_velocity_w=avg_velocity_w,
                  velocity_u_std=velocity_u_std,
                  velocity_v_std=velocity_v_std,
                  velocity_w_std=velocity_w_std)

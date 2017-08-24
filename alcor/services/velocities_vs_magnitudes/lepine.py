from statistics import (mean,
                        stdev)
from typing import (Union,
                    Iterator,
                    Tuple,
                    List)

from alcor.models import (Group,
                          Star)
from alcor.models.velocities_vs_magnitudes import (LepineCaseUCloud,
                                                   LepineCaseVCloud,
                                                   LepineCaseWCloud,
                                                   LepineCaseUBin,
                                                   LepineCaseVBin,
                                                   LepineCaseWBin)
from alcor.types import StarsBinsType
from .utils import (STARS_BIN_SIZE,
                    STARS_BINS_COUNT,
                    DEFAULT_VELOCITY_STD,
                    MIN_BOLOMETRIC_MAGNITUDE,
                    star_bolometric_index)


def clouds(stars: List[Star],
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


def bins(*,
         stars: List[Star],
         group: Group) -> Iterator[Union[LepineCaseUBin,
                                         LepineCaseVBin,
                                         LepineCaseWBin]]:
    u_stars_bins, v_stars_bins, w_stars_bins = stars_bins_by_velocities(stars)
    yield from u_bins(stars_bins=u_stars_bins,
                      group=group)
    yield from v_bins(stars_bins=v_stars_bins,
                      group=group)
    yield from w_bins(stars_bins=w_stars_bins,
                      group=group)


def stars_bins_by_velocities(stars: List[Star]) -> Tuple[StarsBinsType,
                                                         StarsBinsType,
                                                         StarsBinsType]:
    u_stars_bins = [[] for _ in range(STARS_BINS_COUNT)]
    v_stars_bins = [[] for _ in range(STARS_BINS_COUNT)]
    w_stars_bins = [[] for _ in range(STARS_BINS_COUNT)]

    for star in stars:
        index = star_bolometric_index(star)

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


def u_bins(*,
           stars_bins: StarsBinsType,
           group: Group) -> Iterator[LepineCaseUBin]:
    for index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + STARS_BIN_SIZE * (index + 0.5))
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


def v_bins(*,
           stars_bins: StarsBinsType,
           group: Group) -> Iterator[LepineCaseVBin]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + STARS_BIN_SIZE * (stars_bin_index + 0.5))
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


def w_bins(*,
           stars_bins: StarsBinsType,
           group: Group) -> Iterator[LepineCaseWBin]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        if not stars_bin:
            continue

        avg_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                         + STARS_BIN_SIZE * (stars_bin_index + 0.5))
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

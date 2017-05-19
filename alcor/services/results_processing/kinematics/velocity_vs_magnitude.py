import csv
from contextlib import ExitStack
from math import ceil
from statistics import (mean,
                        stdev)
from typing import (List,
                    Tuple)

from collections import Iterable

from alcor.models.star import Star
from alcor.types import (RowType,
                         StarsBinsType)

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)


def write_data_for_velocities_vs_magnitude(stars: List[Star]) -> None:
    # TODO: add reading criterion flag as command line arg
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        write_cloud_data_lepine_case(stars)
        write_bins_data_lepine_case(stars)
    else:
        write_cloud_data_raw_case(stars)
        write_bins_data_raw_case(stars)


def write_cloud_data_lepine_case(stars: List[Star]) -> None:
    with ExitStack() as stack:
        u_vs_mag_cloud_file = stack.enter_context(
            open('u_vs_mag_cloud.csv', mode='w'))
        v_vs_mag_cloud_file = stack.enter_context(
            open('v_vs_mag_cloud.csv', mode='w'))
        w_vs_mag_cloud_file = stack.enter_context(
            open('w_vs_mag_cloud.csv', mode='w'))

        u_vs_mag_cloud_writer = csv.writer(u_vs_mag_cloud_file,
                                           delimiter='  ')
        v_vs_mag_cloud_writer = csv.writer(v_vs_mag_cloud_file,
                                           delimiter='  ')
        w_vs_mag_cloud_writer = csv.writer(w_vs_mag_cloud_file,
                                           delimiter='  ')

        u_vs_mag_cloud_writer.writerow(['bolometric_magnitude',
                                        'velocity_u'])
        v_vs_mag_cloud_writer.writerow(['bolometric_magnitude',
                                        'velocity_v'])
        w_vs_mag_cloud_writer.writerow(['bolometric_magnitude',
                                        'velocity_w'])

        u_vs_mag_stars, v_vs_mag_stars, w_vs_mag_stars = generate_clouds(stars)

        for star in u_vs_mag_stars:
            u_vs_mag_cloud_writer.writerow([star.bolometric_magnitude,
                                            star.velocity_u])
        for star in v_vs_mag_stars:
            v_vs_mag_cloud_writer.writerow([star.bolometric_magnitude,
                                            star.velocity_v])

        for star in w_vs_mag_stars:
            w_vs_mag_cloud_writer.writerow([star.bolometric_magnitude,
                                            star.velocity_w])


def generate_clouds(stars: List[Star]) -> Tuple[List[Star],
                                                List[Star],
                                                List[Star]]:
    u_vs_mag_clouds = []
    v_vs_mag_clouds = []
    w_vs_mag_clouds = []

    for star in stars:
        if star.coordinate_x == highest_coordinate(star):
            u_vs_mag_clouds.append(star)
        elif star.coordinate_y == highest_coordinate(star):
            v_vs_mag_clouds.append(star)
        else:
            w_vs_mag_clouds.append(star)
    return (u_vs_mag_clouds,
            v_vs_mag_clouds,
            w_vs_mag_clouds)


def highest_coordinate(star: Star) -> float:
    return max(star.coordinate_x,
               star.coordinate_y,
               star.coordinate_z)


def write_bins_data_lepine_case(stars: List[Star]) -> None:
    u_vs_mag_bins, v_vs_mag_bins, w_vs_mag_bins = generate_stars_bins(stars)

    # TODO: find the way to shorten the next 3 code blocks
    with open('u_vs_mag_bins.csv', 'w') as u_file:
        u_writer = csv.writer(u_file, delimiter='  ')
        u_writer.writerow(['average_bin_magnitude',
                           'average_velocity_u',
                           'velocity_u_std'])
        for stars_bin_index, stars_bin in enumerate(u_vs_mag_bins):
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (stars_bin_index - 0.5))
            average_bin_velocity_u = mean(star.velocity_u
                                          for star in stars_bin)
            bin_standard_deviation_of_velocity_u = stdev(star.velocity_u
                                                         for star in stars_bin)
            u_writer.writerow([average_bin_magnitude,
                               average_bin_velocity_u,
                               bin_standard_deviation_of_velocity_u])

    with open('v_vs_mag_bins.csv', 'w') as v_file:
        v_writer = csv.writer(v_file, delimiter='  ')
        v_writer.writerow(['average_bin_magnitude',
                           'average_velocity_v',
                           'velocity_v_std'])
        for stars_bin_index, stars_bin in enumerate(v_vs_mag_bins):
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (stars_bin_index - 0.5))
            average_bin_velocity_v = mean(star.velocity_v
                                          for star in stars_bin)
            bin_standard_deviation_of_velocity_v = stdev(star.velocity_v
                                                         for star in stars_bin)
            v_writer.writerow([average_bin_magnitude,
                               average_bin_velocity_v,
                               bin_standard_deviation_of_velocity_v])

    with open('w_vs_mag_bins.csv', 'w') as w_file:
        w_writer = csv.writer(w_file, delimiter='  ')
        w_writer.writerow(['average_bin_magnitude',
                           'average_velocity_w',
                           'velocity_w_std'])
        for stars_bin_index, stars_bin in enumerate(w_vs_mag_bins):
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (stars_bin_index - 0.5))
            average_bin_velocity_w = mean(star.velocity_w
                                          for star in stars_bin)
            bin_standard_deviation_of_velocity_w = stdev(star.velocity_w
                                                         for star in stars_bin)
            w_writer.writerow([average_bin_magnitude,
                               average_bin_velocity_w,
                               bin_standard_deviation_of_velocity_w])


def generate_stars_bins(stars: List[Star]) -> Tuple[StarsBinsType,
                                                    StarsBinsType,
                                                    StarsBinsType]:
    u_vs_mag_bins = [[] for _ in range(BINS_COUNT)]
    v_vs_mag_bins = [[] for _ in range(BINS_COUNT)]
    w_vs_mag_bins = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        stars_bin_index = get_magnitude_bin_index(star)

        if star.coordinate_x == highest_coordinate(star):
            v_vs_mag_bins[stars_bin_index].append(star)
            w_vs_mag_bins[stars_bin_index].append(star)
        elif star.coordinate_y == highest_coordinate(star):
            u_vs_mag_bins[stars_bin_index].append(star)
            w_vs_mag_bins[stars_bin_index].append(star)
        else:
            u_vs_mag_bins[stars_bin_index].append(star)
            v_vs_mag_bins[stars_bin_index].append(star)
    return u_vs_mag_bins, v_vs_mag_bins, w_vs_mag_bins


def get_magnitude_bin_index(star: Star) -> int:
    return int(ceil((star.bolometric_magnitude
                     - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))


def write_cloud_data_raw_case(stars: List[Star]) -> None:
    with open('uvw_vs_mag_cloud.csv', 'w') as uvw_file:
        uvw_writer = csv.writer(uvw_file, delimiter='  ')
        uvw_writer.writerow(['bolometric_magnitude',
                             'velocity_u',
                             'velocity_v',
                             'velocity_w'])
        for star in stars:
            uvw_writer.writerow([star.bolometric_magnitude,
                                 star.velocity_u,
                                 star.velocity_v,
                                 star.velocity_w])


def write_bins_data_raw_case(stars: List[Star]) -> None:
    uvw_vs_mag_bins = [[] for _ in range(BINS_COUNT)]

    for star in stars:
        stars_bin_index = get_magnitude_bin_index(star)
        uvw_vs_mag_bins[stars_bin_index].append(star)

    with open('uvw_vs_mag_bins.csv', 'w') as uvw_file:
        uvw_writer = csv.writer(uvw_file, delimiter='  ')
        uvw_writer.writerow(['average_bin_magnitude',
                             'average_velocity_u',
                             'average_velocity_v',
                             'average_velocity_w',
                             'velocity_u_std',
                             'velocity_v_std',
                             'velocity_w_std'])
        for row in uvw_rows(uvw_vs_mag_bins):
            uvw_writer.writerow(row)


def uvw_rows(stars_bins: StarsBinsType) -> Iterable[RowType]:
    for stars_bin_index, stars_bin in enumerate(stars_bins):
        average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                 + BIN_SIZE * (stars_bin_index
                                               - 0.5))
        average_bin_velocity_u = mean(star.velocity_u
                                      for star in stars_bin)
        average_bin_velocity_v = mean(star.velocity_v
                                      for star in stars_bin)
        average_bin_velocity_w = mean(star.velocity_w
                                      for star in stars_bin)
        bin_standard_deviation_of_velocity_u = stdev(star.velocity_u
                                                     for star in stars_bin)
        bin_standard_deviation_of_velocity_v = stdev(star.velocity_v
                                                     for star in stars_bin)
        bin_standard_deviation_of_velocity_w = stdev(star.velocity_w
                                                     for star in stars_bin)
        yield (average_bin_magnitude,
               average_bin_velocity_u,
               average_bin_velocity_v,
               average_bin_velocity_w,
               bin_standard_deviation_of_velocity_u,
               bin_standard_deviation_of_velocity_v,
               bin_standard_deviation_of_velocity_w)

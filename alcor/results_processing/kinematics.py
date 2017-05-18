import csv
from contextlib import ExitStack
from math import ceil
from statistics import (mean,
                        stdev)
from typing import List

from astropy import units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.coordinates import Distance

from alcor.models.star import Star

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)


def write_data_for_velocity_clouds(stars: List[Star]) -> None:
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        with ExitStack() as stack:
            uv_file = stack.enter_context(open('uv_cloud.csv', mode='w'))
            uw_file = stack.enter_context(open('uw_cloud.csv', mode='w'))
            vw_file = stack.enter_context(open('vw_cloud.csv', mode='w'))

            uv_writer = csv.writer(uv_file, delimiter='  ')
            uw_writer = csv.writer(uw_file, delimiter='  ')
            vw_writer = csv.writer(vw_file, delimiter='  ')

            uv_writer.writerow('velocity_u',
                               'velocity_v')
            uw_writer.writerow('velocity_u',
                               'velocity_w')
            vw_writer.writerow('velocity_v',
                               'velocity_w')

            for star in stars:
                equatorial_coordinates = SkyCoord(
                    ra=star.right_ascension * u.degree,
                    dec=star.declination * u.degree,
                    distance=star.distance * u.kpc)
                coordinate_x = equatorial_coordinates.cartesian.x
                coordinate_y = equatorial_coordinates.cartesian.y
                coordinate_z = equatorial_coordinates.cartesian.z
                highest_coordinate = max(coordinate_x,
                                         coordinate_y,
                                         coordinate_z)
                if coordinate_x == highest_coordinate:
                    vw_writer.writerow(star.velocity_v,
                                       star.velocity_w)
                elif coordinate_y == highest_coordinate:
                    uw_writer.writerow(star.velocity_u,
                                       star.velocity_w)
                elif coordinate_z == highest_coordinate:
                    uv_writer.writerow(star.velocity_u,
                                       star.velocity_v)
            vw_writer.writerow('\n')
            uw_writer.writerow('\n')
            uv_writer.writerow('\n')
    else:
        with open('uvw_cloud.csv', 'w') as uvw_file:
            uvw_writer = csv.writer(uvw_file, delimiter='  ')
            uvw_writer.writerow('velocity_u',
                                'velocity_v',
                                'velocity_w')
            for star in stars:
                uvw_writer.writerow(star.velocity_u,
                                    star.velocity_v,
                                    star.velocity_w)
            uvw_writer.writerow('\n')


def write_data_for_velocities_vs_magnitude(stars: List[Star]) -> None:
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        with ExitStack() as stack:
            u_file = stack.enter_context(open('u_vs_mag_cloud.csv', mode='w'))
            v_file = stack.enter_context(open('v_vs_mag_cloud.csv', mode='w'))
            w_file = stack.enter_context(open('w_vs_mag_cloud.csv', mode='w'))

            u_writer = csv.writer(u_file, delimiter='  ')
            v_writer = csv.writer(v_file, delimiter='  ')
            w_writer = csv.writer(w_file, delimiter='  ')

            u_writer.writerow('bolometric_magnitude',
                              'velocity_u')
            v_writer.writerow('bolometric_magnitude',
                              'velocity_v')
            w_writer.writerow('bolometric_magnitude',
                              'velocity_w')

            u_vs_mag_bins = [[] for _ in range(BINS_COUNT)]
            v_vs_mag_bins = [[] for _ in range(BINS_COUNT)]
            w_vs_mag_bins = [[] for _ in range(BINS_COUNT)]

            for star in stars:
                equatorial_coordinates = SkyCoord(
                    ra=star.right_ascension * u.degree,
                    dec=star.declination * u.degree,
                    distance=star.distance * u.kpc)

                coordinate_x = equatorial_coordinates.cartesian.x
                coordinate_y = equatorial_coordinates.cartesian.y
                coordinate_z = equatorial_coordinates.cartesian.z
                highest_coordinate = max(coordinate_x,
                                         coordinate_y,
                                         coordinate_z)

                stars_bin_index = int(ceil((star.bolometric_magnitude
                                            - MIN_BOLOMETRIC_MAGNITUDE)
                                           / BIN_SIZE))

                if coordinate_x == highest_coordinate:
                    v_writer.writerow(star.bolometric_magnitude,
                                      star.velocity_v)
                    w_writer.writerow(star.bolometric_magnitude,
                                      star.velocity_w)
                    v_vs_mag_bins[stars_bin_index].append(star)
                    w_vs_mag_bins[stars_bin_index].append(star)
                elif coordinate_y == highest_coordinate:
                    u_writer.writerow(star.bolometric_magnitude,
                                      star.velocity_u)
                    w_writer.writerow(star.bolometric_magnitude,
                                      star.velocity_w)
                    u_vs_mag_bins[stars_bin_index].append(star)
                    w_vs_mag_bins[stars_bin_index].append(star)
                elif coordinate_z == highest_coordinate:
                    u_writer.writerow(star.bolometric_magnitude,
                                      star.velocity_u)
                    v_writer.writerow(star.bolometric_magnitude,
                                      star.velocity_v)
                    u_vs_mag_bins[stars_bin_index].append(star)
                    v_vs_mag_bins[stars_bin_index].append(star)

            u_writer.writerow('\n')
            v_writer.writerow('\n')
            w_writer.writerow('\n')

        with open('u_vs_mag_bins.csv', 'w') as u_file:
            u_writer = csv.writer(u_file, delimiter='  ')
            u_writer.writerow('average_bin_magnitude',
                              'average_velocity_u',
                              'velocity_u_std')
            for stars_bin_index, stars_bin in enumerate(u_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (stars_bin_index - 0.5))
                average_bin_velocity_u = mean(star.velocity_u
                                              for star in stars_bin)
                bin_standard_deviation_of_velocity_u = stdev(
                    star.velocity_u
                    for star in stars_bin)
                u_writer.writerow(average_bin_magnitude,
                                  average_bin_velocity_u,
                                  bin_standard_deviation_of_velocity_u)
            u_writer.writerow('\n')

        with open('v_vs_mag_bins.csv', 'w') as v_file:
            v_writer = csv.writer(v_file, delimiter='  ')
            v_writer.writerow('average_bin_magnitude',
                              'average_velocity_v',
                              'velocity_v_std')
            for stars_bin_index, stars_bin in enumerate(v_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (stars_bin_index - 0.5))
                average_bin_velocity_v = mean(star.velocity_v
                                              for star in stars_bin)
                bin_standard_deviation_of_velocity_v = stdev(
                    star.velocity_v
                    for star in stars_bin)
                v_writer.writerow(average_bin_magnitude,
                                  average_bin_velocity_v,
                                  bin_standard_deviation_of_velocity_v)
            v_writer.writerow('\n')

        with open('w_vs_mag_bins.csv', 'w') as w_file:
            w_writer = csv.writer(w_file, delimiter='  ')
            w_writer.writerow('average_bin_magnitude',
                              'average_velocity_w',
                              'velocity_w_std')
            for stars_bin_index, stars_bin in enumerate(w_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (stars_bin_index - 0.5))
                average_bin_velocity_w = mean(star.velocity_w
                                              for star in stars_bin)
                bin_standard_deviation_of_velocity_w = stdev(
                    star.velocity_w
                    for star in stars_bin)
                w_writer.writerow(average_bin_magnitude,
                                  average_bin_velocity_w,
                                  bin_standard_deviation_of_velocity_w)
            w_writer.writerow('\n')
    else:
        with open('uvw_vs_mag_cloud.csv', 'w') as uvw_file:
            uvw_writer = csv.writer(uvw_file, delimiter='  ')
            uvw_writer.writerow('bolometric_magnitude',
                                'velocity_u',
                                'velocity_v',
                                'velocity_w')

            uvw_vs_mag_bins = [[] for _ in range(BINS_COUNT)]

            for star in stars:
                stars_bin_index = int(ceil((star.bolometric_magnitude
                                            - MIN_BOLOMETRIC_MAGNITUDE)
                                           / BIN_SIZE))
                uvw_writer.writerow(star.bolometric_magnitude,
                                    star.velocity_u,
                                    star.velocity_v,
                                    star.velocity_w)
                uvw_vs_mag_bins[stars_bin_index].append(star)
            uvw_writer.writerow('\n')

        with open('uvw_vs_mag_bins.csv', 'w') as uvw_file:
            uvw_writer = csv.writer(uvw_file, delimiter='  ')
            uvw_writer.writerow('average_bin_magnitude',
                                'average_velocity_u',
                                'average_velocity_v',
                                'average_velocity_w',
                                'velocity_u_std',
                                'velocity_v_std',
                                'velocity_w_std')
            for stars_bin_index, stars_bin in enumerate(uvw_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (stars_bin_index
                                                       - 0.5))
                average_bin_velocity_u = mean(star.velocity_u
                                              for star in stars_bin)
                average_bin_velocity_v = mean(star.velocity_v
                                              for star in stars_bin)
                average_bin_velocity_w = mean(star.velocity_w
                                              for star in stars_bin)
                bin_standard_deviation_of_velocity_u = stdev(
                    star.velocity_u
                    for star in stars_bin)
                bin_standard_deviation_of_velocity_v = stdev(
                    star.velocity_v
                    for star in stars_bin)
                bin_standard_deviation_of_velocity_w = stdev(
                    star.velocity_w
                    for star in stars_bin)
                uvw_writer.writerow(average_bin_magnitude,
                                    average_bin_velocity_u,
                                    average_bin_velocity_v,
                                    average_bin_velocity_w,
                                    bin_standard_deviation_of_velocity_u,
                                    bin_standard_deviation_of_velocity_v,
                                    bin_standard_deviation_of_velocity_w)
            uvw_writer.writerow('\n')

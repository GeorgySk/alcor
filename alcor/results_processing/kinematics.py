import csv
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
        with open('uv_cloud.csv', 'w') as uv_file:
            uv_writer = csv.writer(uv_file, delimiter='  ')
            uv_writer.writerow('star.velocity_u',
                               'star.velocity_v')
            with open('uw_cloud.csv', 'w') as uw_file:
                uw_writer = csv.writer(uw_file, delimiter='  ')
                uw_writer.writerow('star.velocity_u',
                                   'star.velocity_w')
                with open('vw_cloud.csv', 'w') as vw_file:
                    vw_writer = csv.writer(vw_file, delimiter='  ')
                    vw_writer.writerow('star.velocity_v',
                                       'star.velocity_w')

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
                        else:
                            print('Game over')
                    vw_writer.writerow('\n')
                    uw_writer.writerow('\n')
                    uv_writer.writerow('\n')
    else:
        with open('uvw_cloud.csv', 'w') as uvw_file:
            uvw_writer = csv.writer(uvw_file, delimiter='  ')
            uvw_writer.writerow('star.velocity_u',
                                'star.velocity_v',
                                'star.velocity_w')
            for star in stars:
                uvw_writer.writerow(star.velocity_u,
                                    star.velocity_v,
                                    star.velocity_w)
            uvw_writer.writerow('\n')


def write_data_for_velocities_vs_magnitude(stars: List[Star]) -> None:
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        with open('u_vs_mag_cloud.csv', 'w') as u_file:
            u_writer = csv.writer(u_file, delimiter='  ')
            u_writer.writerow('bolometric_magnitude',
                              'star.velocity_u')
            with open('v_vs_mag_cloud.csv', 'w') as v_file:
                v_writer = csv.writer(v_file, delimiter='  ')
                v_writer.writerow('bolometric_magnitude',
                                  'star.velocity_v')
                with open('w_vs_mag_cloud.csv', 'w') as w_file:
                    w_writer = csv.writer(w_file, delimiter='  ')
                    w_writer.writerow('bolometric_magnitude',
                                      'star.velocity_w')

                    u_vs_mag_bins = []
                    v_vs_mag_bins = []
                    w_vs_mag_bins = []

                    for star in stars:
                        bolometric_magnitude = 2.5 * star.luminosity + 4.75

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
                        data_bin = int(ceil((bolometric_magnitude
                                             - MIN_BOLOMETRIC_MAGNITUDE)
                                            / BIN_SIZE))
                        if coordinate_x == highest_coordinate:
                            v_writer.writerow(bolometric_magnitude,
                                              star.velocity_v)
                            w_writer.writerow(bolometric_magnitude,
                                              star.velocity_w)
                            v_vs_mag_bins[data_bin].append(star)
                            w_vs_mag_bins[data_bin].append(star)
                        elif coordinate_y == highest_coordinate:
                            u_writer.writerow(bolometric_magnitude,
                                              star.velocity_u)
                            w_writer.writerow(bolometric_magnitude,
                                              star.velocity_w)
                            u_vs_mag_bins[data_bin].append(star)
                            w_vs_mag_bins[data_bin].append(star)
                        elif coordinate_z == highest_coordinate:
                            u_writer.writerow(bolometric_magnitude,
                                              star.velocity_u)
                            v_writer.writerow(bolometric_magnitude,
                                              star.velocity_v)
                            u_vs_mag_bins[data_bin].append(star)
                            v_vs_mag_bins[data_bin].append(star)
                        else:
                            print('Game over')
                    u_writer.writerow('\n')
                    v_writer.writerow('\n')
                    w_writer.writerow('\n')

        with open('u_vs_mag_bins.csv', 'w') as u_file:
            u_writer = csv.writer(u_file, delimiter='  ')
            u_writer.writerow('average_bin_magnitude',
                              'average_bin_velocity_u',
                              'bin_standard_deviation_of_velocity_u')
            for data_bin_index, data_bin in enumerate(u_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (data_bin_index - 0.5))
                average_bin_velocity_u = mean(star.velocity_u
                                              for star in data_bin)
                bin_standard_deviation_of_velocity_u = stdev(star.velocity_u
                                                             for star in
                                                             data_bin)
                u_writer.writerow(average_bin_magnitude,
                                  average_bin_velocity_u,
                                  bin_standard_deviation_of_velocity_u)
            u_writer.writerow('\n')

        with open('v_vs_mag_bins.csv', 'w') as v_file:
            v_writer = csv.writer(v_file, delimiter='  ')
            v_writer.writerow('average_bin_magnitude',
                              'average_bin_velocity_v',
                              'bin_standard_deviation_of_velocity_v')
            for data_bin_index, data_bin in enumerate(v_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (data_bin_index - 0.5))
                average_bin_velocity_v = mean(star.velocity_v
                                              for star in data_bin)
                bin_standard_deviation_of_velocity_v = stdev(star.velocity_v
                                                             for star in
                                                             data_bin)
                v_writer.writerow(average_bin_magnitude,
                                  average_bin_velocity_v,
                                  bin_standard_deviation_of_velocity_v)
            v_writer.writerow('\n')

        with open('w_vs_mag_bins.csv', 'w') as w_file:
            w_writer = csv.writer(w_file, delimiter='  ')
            w_writer.writerow('average_bin_magnitude',
                              'average_bin_velocity_w',
                              'bin_standard_deviation_of_velocity_w')
            for data_bin_index, data_bin in enumerate(w_vs_mag_bins):
                average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                         + BIN_SIZE * (data_bin_index - 0.5))
                average_bin_velocity_w = mean(star.velocity_w
                                              for star in data_bin)
                bin_standard_deviation_of_velocity_w = stdev(star.velocity_w
                                                             for star in
                                                             data_bin)
                w_writer.writerow(average_bin_magnitude,
                                  average_bin_velocity_w,
                                  bin_standard_deviation_of_velocity_w)
            w_writer.writerow('\n')
    else:
        with open('uvw_vs_mag_cloud.csv', 'w') as uvw_file:
            uvw_writer = csv.writer(uvw_file, delimiter='  ')
            uvw_writer.writerow('bolometric_magnitude',
                                'star.velocity_u',
                                'star.velocity_v',
                                'star.velocity_w')

            uvw_vs_mag_bins = []

            for star in stars:
                bolometric_magnitude = 2.5 * star.luminosity + 4.75
                data_bin = int(ceil((bolometric_magnitude
                                     - MIN_BOLOMETRIC_MAGNITUDE)
                                    / BIN_SIZE))
                uvw_writer.writerow(bolometric_magnitude,
                                    star.velocity_u,
                                    star.velocity_v,
                                    star.velocity_w)
                uvw_vs_mag_bins[data_bin].append(star)
            uvw_writer.writerow('\n')

            with open('uvw_vs_mag_bins.csv.csv', 'w') as uvw_file:
                uvw_writer = csv.writer(uvw_file, delimiter='  ')
                uvw_writer.writerow('average_bin_magnitude',
                                    'average_bin_velocity_u',
                                    'average_bin_velocity_v',
                                    'average_bin_velocity_w',
                                    'bin_standard_deviation_of_velocity_u',
                                    'bin_standard_deviation_of_velocity_v',
                                    'bin_standard_deviation_of_velocity_w')
                for data_bin_index, data_bin in enumerate(uvw_vs_mag_bins):
                    average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                             + BIN_SIZE * (data_bin_index
                                                           - 0.5))
                    average_bin_velocity_u = mean(
                        star.velocity_u for star in data_bin)
                    average_bin_velocity_v = mean(
                        star.velocity_v for star in data_bin)
                    average_bin_velocity_w = mean(
                        star.velocity_w for star in data_bin)
                    bin_standard_deviation_of_velocity_u = stdev(
                        star.velocity_u
                        for star in data_bin)
                    bin_standard_deviation_of_velocity_v = stdev(
                        star.velocity_v
                        for star in data_bin)
                    bin_standard_deviation_of_velocity_w = stdev(
                        star.velocity_w
                        for star in data_bin)
                    uvw_writer.writerow(average_bin_magnitude,
                                        average_bin_velocity_u,
                                        average_bin_velocity_v,
                                        average_bin_velocity_w,
                                        bin_standard_deviation_of_velocity_u,
                                        bin_standard_deviation_of_velocity_v,
                                        bin_standard_deviation_of_velocity_w)
                uvw_writer.writerow('\n')

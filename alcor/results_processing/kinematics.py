import csv
from statistics import mean, stdev
from typing import List

from astropy import units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.coordinates import Distance

from alcor.models.star import Star
from alcor.results_processing.luminosity_function import (
    MIN_BOLOMETRIC_MAGNITUDE,
    BIN_SIZE)


def write_bins_kinematic_info(data_bins: List[List[Star]]) -> None:
    with open('magnitude_bins.csv', 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='  ')
        for data_bin_index, data_bin in enumerate(data_bins):
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (data_bin_index - 0.5))
            average_bin_velocity_u = mean(star.velocity_u for star in data_bin)
            average_bin_velocity_v = mean(star.velocity_v for star in data_bin)
            average_bin_velocity_w = mean(star.velocity_w for star in data_bin)
            bin_standard_deviation_of_velocity_u = stdev(star.velocity_u
                                                         for star in data_bin)
            bin_standard_deviation_of_velocity_v = stdev(star.velocity_v
                                                         for star in data_bin)
            bin_standard_deviation_of_velocity_w = stdev(star.velocity_w
                                                         for star in data_bin)
            output_writer.writerow(average_bin_magnitude,
                                   average_bin_velocity_u,
                                   average_bin_velocity_v,
                                   average_bin_velocity_w,
                                   bin_standard_deviation_of_velocity_u,
                                   bin_standard_deviation_of_velocity_v,
                                   bin_standard_deviation_of_velocity_w)
        output_writer.writerow('\n')


def write_data_for_velocity_clouds(stars: List[Star]) -> None:
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        with open('uv_cloud.csv', 'w') as uv_file:
            uv_writer = csv.writer(uv_file, delimiter='  ')
            with open('uw_cloud.csv', 'w') as uw_file:
                uw_writer = csv.writer(uw_file, delimiter='  ')
                with open('vw_cloud.csv', 'w') as vw_file:
                    vw_writer = csv.writer(vw_file, delimiter='  ')

                    for star in stars:
                        equatorial_coordinates \
                            = SkyCoord(ra=star.right_ascension * u.degree,
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
            with open('v_vs_mag_cloud.csv', 'w') as v_file:
                v_writer = csv.writer(v_file, delimiter='  ')
                with open('w_vs_mag_cloud.csv', 'w') as w_file:
                    w_writer = csv.writer(w_file, delimiter='  ')

                    for star in stars:
                        bolometric_magnitude = 2.5 * star.luminosity + 4.75

                        equatorial_coordinates \
                            = SkyCoord(ra=star.right_ascension * u.degree,
                                       dec=star.declination * u.degree,
                                       distance=star.distance * u.kpc)
                        coordinate_x = equatorial_coordinates.cartesian.x
                        coordinate_y = equatorial_coordinates.cartesian.y
                        coordinate_z = equatorial_coordinates.cartesian.z
                        highest_coordinate = max(coordinate_x,
                                                 coordinate_y,
                                                 coordinate_z)
                        if coordinate_x == highest_coordinate:
                            v_writer.writerow(bolometric_magnitude,
                                              star.velocity_v)
                            w_writer.writerow(bolometric_magnitude,
                                              star.velocity_w)
                        elif coordinate_y == highest_coordinate:
                            u_writer.writerow(bolometric_magnitude,
                                              star.velocity_u)
                            w_writer.writerow(bolometric_magnitude,
                                              star.velocity_w)
                        elif coordinate_z == highest_coordinate:
                            u_writer.writerow(bolometric_magnitude,
                                              star.velocity_u)
                            v_writer.writerow(bolometric_magnitude,
                                              star.velocity_v)
                        else:
                            print('Game over')
                    u_writer.writerow('\n')
                    v_writer.writerow('\n')
                    w_writer.writerow('\n')
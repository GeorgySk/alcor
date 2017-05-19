import csv
from contextlib import ExitStack
from typing import List

from astropy import units as u
from astropy.coordinates.sky_coordinate import SkyCoord

from alcor.models.star import Star


def write_data_for_velocity_clouds(stars: List[Star]) -> None:
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        with ExitStack() as stack:
            uv_cloud_file = stack.enter_context(open('uv_cloud.csv', mode='w'))
            uw_cloud_file = stack.enter_context(open('uw_cloud.csv', mode='w'))
            vw_cloud_file = stack.enter_context(open('vw_cloud.csv', mode='w'))

            uv_cloud_writer = csv.writer(uv_cloud_file, delimiter='  ')
            uw_cloud_writer = csv.writer(uw_cloud_file, delimiter='  ')
            vw_cloud_writer = csv.writer(vw_cloud_file, delimiter='  ')

            uv_cloud_writer.writerow(['velocity_u',
                                      'velocity_v'])
            uw_cloud_writer.writerow(['velocity_u',
                                      'velocity_w'])
            vw_cloud_writer.writerow(['velocity_v',
                                      'velocity_w'])

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
                    vw_cloud_writer.writerow([star.velocity_v,
                                              star.velocity_w])
                elif coordinate_y == highest_coordinate:
                    uw_cloud_writer.writerow([star.velocity_u,
                                              star.velocity_w])
                elif coordinate_z == highest_coordinate:
                    uv_cloud_writer.writerow([star.velocity_u,
                                              star.velocity_v])
    else:
        with open('uvw_cloud.csv', 'w') as uvw_file:
            uvw_writer = csv.writer(uvw_file, delimiter='  ')
            uvw_writer.writerow(['velocity_u',
                                 'velocity_v',
                                 'velocity_w'])
            for star in stars:
                uvw_writer.writerow([star.velocity_u,
                                     star.velocity_v,
                                     star.velocity_w])

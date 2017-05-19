import csv
from contextlib import ExitStack
from typing import List

from astropy import units as u
from astropy.coordinates.sky_coordinate import SkyCoord

from alcor.models.star import Star


def write_data_for_velocity_clouds(stars: List[Star]) -> None:
    lepine_selection_criterion_applied = True

    if lepine_selection_criterion_applied:
        write_data_lepin_case(stars)
    else:
        write_data_raw_case(stars)


def write_data_lepin_case(stars: List[Star]) -> None:
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
            if star.coordinate_x == highest_coordinate(star):
                vw_cloud_writer.writerow([star.velocity_v,
                                          star.velocity_w])
            elif star.coordinate_y == highest_coordinate(star):
                uw_cloud_writer.writerow([star.velocity_u,
                                          star.velocity_w])
            else:
                uv_cloud_writer.writerow([star.velocity_u,
                                          star.velocity_v])


# TODO: this function is also in velocity_vs_magnitude. Where should I put it?
def highest_coordinate(star: Star) -> float:
    return max(star.coordinate_x,
               star.coordinate_y,
               star.coordinate_z)


def write_data_raw_case(stars: List[Star]) -> None:
    with open('uvw_cloud.csv', 'w') as uvw_file:
        uvw_writer = csv.writer(uvw_file, delimiter='  ')
        uvw_writer.writerow(['velocity_u',
                             'velocity_v',
                             'velocity_w'])
        for star in stars:
            uvw_writer.writerow([star.velocity_u,
                                 star.velocity_v,
                                 star.velocity_w])

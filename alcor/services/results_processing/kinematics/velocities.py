import csv
from contextlib import ExitStack
from typing import (List,
                    Tuple)

from alcor.models.star import Star
from alcor.models.velocities import LepineCaseUVCloud


def write_velocity_clouds_data(stars: List[Star],
                               lepine_criterion_applied: bool) -> None:
    if lepine_criterion_applied:
        write_data_lepin_case(stars)
    else:
        write_data_raw_case(stars)


def write_data_lepin_case(stars: List[Star]) -> None:
    uv_cloud_stars, uw_cloud_stars, vw_cloud_stars = generate_clouds(stars)

    with ExitStack() as stack:
        uv_cloud_file = stack.enter_context(open('uv_cloud.csv', mode='w'))
        uw_cloud_file = stack.enter_context(open('uw_cloud.csv', mode='w'))
        vw_cloud_file = stack.enter_context(open('vw_cloud.csv', mode='w'))

        uv_cloud_writer = csv.writer(uv_cloud_file, delimiter=' ')
        uw_cloud_writer = csv.writer(uw_cloud_file, delimiter=' ')
        vw_cloud_writer = csv.writer(vw_cloud_file, delimiter=' ')

        uv_cloud_writer.writerow(['velocity_u',
                                  'velocity_v'])
        uw_cloud_writer.writerow(['velocity_u',
                                  'velocity_w'])
        vw_cloud_writer.writerow(['velocity_v',
                                  'velocity_w'])

        for star in uv_cloud_stars:
            uv_cloud_writer.writerow([star.velocity_u,
                                      star.velocity_v])
        for star in uw_cloud_stars:
            uw_cloud_writer.writerow([star.velocity_u,
                                      star.velocity_w])
        for star in vw_cloud_stars:
            vw_cloud_writer.writerow([star.velocity_v,
                                      star.velocity_w])


def generate_clouds(stars: List[Star]) -> Tuple[List[Star],
                                                List[Star],
                                                List[Star]]:
    uv_cloud = []
    uw_cloud = []
    vw_cloud = []
    for star in stars:
        highest_coordinate = max(star.coordinate_x,
                                 star.coordinate_y,
                                 star.coordinate_z)
        if star.coordinate_x == highest_coordinate:
            vw_cloud.append(star)
        elif star.coordinate_y == highest_coordinate:
            uw_cloud.append(star)
        else:
            uv_cloud.append(star)
    return uv_cloud, uw_cloud, vw_cloud


def write_data_raw_case(stars: List[Star]) -> None:
    with open('uvw_cloud.csv', 'w') as uvw_file:
        uvw_writer = csv.writer(uvw_file, delimiter=' ')
        uvw_writer.writerow(['velocity_u',
                             'velocity_v',
                             'velocity_w'])
        for star in stars:
            uvw_writer.writerow([star.velocity_u,
                                 star.velocity_v,
                                 star.velocity_w])

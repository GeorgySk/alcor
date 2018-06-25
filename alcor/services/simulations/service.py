import logging
import os
from subprocess import check_call
from typing import Dict

from alcor.services.common import output_file_name
from alcor.types import (GridParametersInfoType,
                         CSVParametersInfoType)
from . import grid

logger = logging.getLogger(__name__)


def run(*,
        geometry: str,
        precision: int,
        grid_parameters_info: GridParametersInfoType,
        csv_parameters_info: CSVParametersInfoType) -> None:
    for parameters_values in grid.parameters_values(
            parameters_info=grid_parameters_info,
            precision=precision):
        run_simulation(parameters_values=parameters_values,
                       csv_parameters_info=csv_parameters_info,
                       geometry=geometry,
                       output_file_name=output_file_name())


def run_simulation(*,
                   parameters_values: Dict[str, float],
                   csv_parameters_info: CSVParametersInfoType,
                   geometry: str,
                   output_file_name: str) -> None:
    args = ['./main.e',
            '-db', parameters_values['DB_fraction'],
            '-g', parameters_values['thin_disk_age'],
            '-tda', parameters_values['thick_disk_age'],
            '-ha', parameters_values['halo_age'],
            '-hsft', parameters_values['halo_stars_formation_time'],
            '-tde', parameters_values['thick_disk_star_formation_exponent'],
            '-mf', parameters_values['initial_mass_function_exponent'],
            '-ifr', parameters_values['lifetime_mass_ratio'],
            '-bt', parameters_values['burst_time'],
            '-mr', parameters_values['mass_reduction_factor'],
            '-tdsf', parameters_values['thick_disk_stars_fraction'],
            '-hsf', parameters_values['halo_stars_fraction'],
            '-rad', parameters_values['radius'],
            '-o', output_file_name,
            '-geom', geometry]

    if geometry == 'cones':
        try:
            args.extend(['-cl', parameters_values['longitudes']])
        except KeyError:
            longitudes_info = csv_parameters_info['longitudes']
            args.extend(['-clcsv', os.path.abspath(longitudes_info['path']),
                         '-clcol', longitudes_info['column']])

        try:
            args.extend(['-cb', parameters_values['latitudes']])
        except KeyError:
            latitudes_dict = csv_parameters_info['latitudes']
            args.extend(['-cbcsv', os.path.abspath(latitudes_dict['path']),
                         '-cbcol', latitudes_dict['column']])

    args = list(map(str, args))
    args_str = ' '.join(args)
    logger.info(f'Invoking simulation with command "{args_str}".')
    check_call(args)

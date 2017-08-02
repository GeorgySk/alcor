import logging
import os
import uuid
from decimal import Decimal
from subprocess import check_call
from typing import (Any,
                    Iterable,
                    Dict)

from sqlalchemy.orm.session import Session

from alcor.models import Group
from alcor.models.simulation import Parameter
from alcor.services.parameters import generate_parameters_values
from alcor.services.restrictions import (OUTPUT_FILE_EXTENSION,
                                         MAX_OUTPUT_FILE_NAME_LENGTH)
from alcor.types import NumericType
from alcor.utils import parse_stars

logger = logging.getLogger(__name__)


def run_simulations(*,
                    settings: Dict[str, Any],
                    session: Session) -> None:
    model_type = settings['model_type']
    geometry = settings['geometry']
    precision = settings['precision']
    parameters_info = settings['parameters']

    for parameters_values in generate_parameters_values(
            parameters=parameters_info,
            precision=precision,
            geometry=geometry):
        group_id = uuid.uuid4()
        group = Group(id=group_id,
                      original_unprocessed_group_id=None)

        parameters = generate_parameters(values=parameters_values,
                                         group=group)

        output_file_name = generate_output_file_name(group_id=str(group_id))

        run_simulation(parameters_values=parameters_values,
                       model_type=model_type,
                       geometry=geometry,
                       output_file_name=output_file_name)

        with open(output_file_name) as output_file:
            stars = list(parse_stars(output_file,
                                     group=group))
        os.remove(output_file_name)

        session.add(group)
        session.add_all(parameters)
        session.add_all(stars)
        session.commit()

        # TODO: all cones stars for 1 experiment should be recorded to 1 group
        if geometry == 'cone':
            with open(file='../test_project/processed_cones.txt',
                      mode='a') as file:
                for parameter in parameters:
                    if parameter.name == 'longitude':
                        longitude = parameter.value
                    if parameter.name == 'latitude':
                        latitude = parameter.value
                row = str(longitude) + ' ' + str(latitude) + '\n'
                file.write(row)


def generate_parameters(*,
                        values: Dict[str, Decimal],
                        group: Group) -> Iterable[Parameter]:
    for parameter_name, parameter_value in values.items():
        yield Parameter(group_id=group.id,
                        name=parameter_name,
                        value=parameter_value)


def run_simulation(*,
                   parameters_values: Dict[str, NumericType],
                   model_type: int,
                   geometry: str,
                   output_file_name: str) -> None:
    args = ['./main.e',
            '-db', parameters_values['DB_fraction'],
            '-g', parameters_values['galaxy_age'],
            '-mf', parameters_values['initial_mass_function_exponent'],
            '-ifr', parameters_values['lifetime_mass_ratio'],
            '-bt', parameters_values['burst_time'],
            '-mr', parameters_values['mass_reduction_factor'],
            '-km', model_type,
            '-o', output_file_name,
            '-geom', geometry]
    if geometry == 'cone':
        args.extend(['-cl', parameters_values['longitude'],
                     '-cb', parameters_values['latitude']])
    args = list(map(str, args))
    args_str = ' '.join(args)
    logger.info(f'Invoking simulation with command "{args_str}".')
    check_call(args)


def generate_output_file_name(group_id: str) -> str:
    base_name = group_id[:MAX_OUTPUT_FILE_NAME_LENGTH]
    return ''.join([base_name, OUTPUT_FILE_EXTENSION])

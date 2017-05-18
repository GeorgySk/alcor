import uuid
import logging

from decimal import Decimal
from subprocess import check_call
from typing import (Dict,
                    Any)

from cassandra.cluster.Session import Session

from alcor.models.parameter import Parameter
from alcor.services.data_access import insert
from alcor.services.parameters import generate_parameters_values
from alcor.types import NumericType
from alcor.utils import parse_stars

OUTPUT_FILE_EXTENSION = '.res'
MAX_OUTPUT_FILE_NAME_LENGTH = 5

logger = logging.getLogger(__name__)


def run_simulations(*,
                    settings: Dict[str, Any],
                    session: Session) -> None:
    model_type = settings['model_type']
    precision = settings['precision']
    parameters_info = settings['parameters']
    for parameters_values in generate_parameters_values(
            parameters_info=parameters_info,
            precision=precision):
        parameters_group_id = uuid.uuid4()
        save_parameters(values=parameters_values,
                        group_id=parameters_group_id,
                        session=session)

        output_file_name = generate_output_file_name(
            parameters_group_id=str(parameters_group_id))

        run_simulation(parameters_values=parameters_values,
                       model_type=model_type,
                       output_file_name=output_file_name)

        save_stars(file_name=output_file_name,
                   group_id=parameters_group_id,
                   session=session)


def save_parameters(*,
                    values: Dict[str, Decimal],
                    group_id: uuid.UUID,
                    session: Session) -> None:
    instances = [Parameter(group_id=group_id,
                           name=parameter_name,
                           value=parameter_value)
                 for parameter_name, parameter_value in values.items()]
    insert(instances=instances,
           session=session)


def run_simulation(*,
                   parameters_values: Dict[str, NumericType],
                   model_type: int,
                   output_file_name: str) -> None:
    args = ['./main.e',
            '-db', parameters_values['DB_fraction'],
            '-g', parameters_values['galaxy_age'],
            '-mf', parameters_values['initial_mass_function_exponent'],
            '-ifr', parameters_values['lifetime_mass_ratio'],
            '-bt', parameters_values['burst_time'],
            '-mr', parameters_values['mass_reduction_factor'],
            '-km', model_type,
            '-o', output_file_name]
    args = list(map(str, args))
    args_str = ' '.join(args)
    logger.info(f'Invoking simulation with command "{args_str}".')
    check_call(args)


def save_stars(file_name: str,
               group_id: uuid.UUID,
               session: Session) -> None:
    with open(file_name) as output_file:
        stars = parse_stars(output_file,
                            group_id=group_id)
        insert(instances=stars,
               session=session)


def generate_output_file_name(parameters_group_id: str
                              ) -> str:
    base_name = parameters_group_id[:MAX_OUTPUT_FILE_NAME_LENGTH]
    return ''.join([base_name, OUTPUT_FILE_EXTENSION])

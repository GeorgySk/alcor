import logging
import uuid
from decimal import Decimal
from subprocess import check_call
from typing import (Any,
                    Iterable,
                    Dict,
                    Tuple)

from cassandra.cluster import Session as CassandraSession
from sqlalchemy.orm.session import Session

from alcor.cassandra_models import CGroup, CStar
from alcor.cassandra_models.simulation import CParameter
from alcor.models import Group
from alcor.models.simulation import Parameter
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from alcor.services.parameters import generate_parameters_values
from alcor.services.restrictions import (OUTPUT_FILE_EXTENSION,
                                         MAX_OUTPUT_FILE_NAME_LENGTH)
from alcor.types import NumericType
from alcor.utils import parse_stars

logger = logging.getLogger(__name__)


def run_simulations(*,
                    settings: Dict[str, Any],
                    cassandra_session: CassandraSession,
                    session: Session) -> None:
    model_type = settings['model_type']
    precision = settings['precision']
    parameters_info = settings['parameters']
    for parameters_values in generate_parameters_values(
            parameters_info=parameters_info,
            precision=precision):
        group_id = uuid.uuid4()
        c_group = CGroup(id=group_id)
        group = Group(id=group_id)

        c_parameters, parameters = zip(*generate_parameters(
            values=parameters_values,
            c_group=c_group,
            group=group))

        output_file_name = generate_output_file_name(group_id=str(group_id))

        run_simulation(parameters_values=parameters_values,
                       model_type=model_type,
                       output_file_name=output_file_name)

        with open(output_file_name) as output_file:
            c_stars, stars = zip(*parse_stars(output_file,
                                              group=group,
                                              c_group=c_group))

        insert_groups_statement = model_insert_statement(CGroup)
        insert(instances=[c_group],
               statement=insert_groups_statement,
               session=cassandra_session)

        insert_parameters_statement = model_insert_statement(CParameter)
        insert(instances=c_parameters,
               statement=insert_parameters_statement,
               session=cassandra_session)

        insert_stars_statement = model_insert_statement(CStar)
        insert(instances=c_stars,
               statement=insert_stars_statement,
               session=cassandra_session)

        session.add(group)
        session.add_all(parameters)
        session.add_all(stars)
        session.commit()


def generate_parameters(*,
                        values: Dict[str, Decimal],
                        c_group: CGroup,
                        group: Group) -> Iterable[Tuple[CParameter,
                                                        Parameter]]:
    for parameter_name, parameter_value in values.items():
        yield (CParameter(group_id=c_group.id,
                          name=parameter_name,
                          value=parameter_value),
               Parameter(group_id=group.id,
                         name=parameter_name,
                         value=parameter_value))


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


def generate_output_file_name(group_id: str) -> str:
    base_name = group_id[:MAX_OUTPUT_FILE_NAME_LENGTH]
    return ''.join([base_name, OUTPUT_FILE_EXTENSION])

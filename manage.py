#!/usr/bin/env python3
import logging
import os
import uuid
from collections import OrderedDict
from decimal import Decimal
from subprocess import check_call
from typing import (Any,
                    Iterable,
                    Dict)

import click
from cassandra.cluster import (Cluster,
                               Session)
from cassandra.cqlengine import connection
from cassandra_helpers.connectable import check_connection
from cassandra_helpers.keyspace import (keyspace_exists,
                                        create_keyspace,
                                        drop_keyspace)
from cassandra_helpers.models import sync_tables

from alcor.models import (STAR_PARAMETERS_NAMES,
                          Parameter,
                          Star)
from alcor.services.data_access import insert
from alcor.services.parameters import generate_parameters_values
from alcor.types import NumericType
from alcor.utils import load_settings

OUTPUT_FILE_EXTENSION = '.res'
MAX_OUTPUT_FILE_NAME_LENGTH = 5

logger = logging.getLogger(__name__)


@click.group()
@click.pass_context
def main(ctx: click.Context) -> None:
    logging.basicConfig(level=logging.DEBUG)
    contact_points = [os.environ['CASSANDRA_RPC_ADDRESS']]
    port = int(os.environ['CASSANDRA_RPC_PORT'])
    ctx.obj = dict(contact_points=contact_points,
                   port=port)


@main.command()
@click.option('--settings-path', '-p',
              default='settings.yml',
              type=click.Path(),
              help='Settings file path '
                   '(absolute or relative, '
                   'default "settings.yml").')
@click.option('--project-dir', '-d',
              required=True,
              type=click.Path(),
              help='Project directory path '
                   '(absolute or relative).')
@click.option('--clean',
              is_flag=True,
              help='Cleans destination keyspace.')
@click.pass_context
def run(ctx: click.Context,
        settings_path: str,
        project_dir: str,
        clean: bool) -> None:
    cluster_settings = ctx.obj
    contact_points = cluster_settings['contact_points']
    port = cluster_settings['port']

    keyspace_name = 'test'

    check_connection(contact_points=contact_points,
                     port=port)
    with Cluster(contact_points=contact_points,
                 port=port) as cluster:
        session = cluster.connect()

        if clean:
            clean_keyspace(keyspace_name,
                           session=session)

        init_db(keyspace_name=keyspace_name,
                session=session)

        settings = load_settings(settings_path)
        os.chdir(project_dir)
        run_simulations(settings=settings,
                        session=session)


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

        simulate(parameters_values=parameters_values,
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


def simulate(*,
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


def parse_stars(lines: Iterable[str],
                group_id: uuid.UUID
                ) -> Iterable[Star]:
    for line in lines:
        parts = line.split()
        params = map(Decimal, parts)
        values = OrderedDict(zip(STAR_PARAMETERS_NAMES,
                                 params))
        yield Star(group_id=group_id,
                   **values)


def generate_output_file_name(parameters_group_id: str
                              ) -> str:
    base_name = parameters_group_id[:MAX_OUTPUT_FILE_NAME_LENGTH]
    return ''.join([base_name, OUTPUT_FILE_EXTENSION])


def init_db(*,
            keyspace_name: str,
            session: Session) -> None:
    init_keyspace(keyspace_name,
                  session=session)
    session.set_keyspace(keyspace_name)
    connection.set_session(session)
    sync_tables(Parameter, Star)


def init_keyspace(name: str,
                  *,
                  session: Session) -> None:
    logger.info(f'Creating "{name}" keyspace.')
    create_keyspace(name,
                    session=session,
                    check_first=True)


def clean_keyspace(name: str,
                   *,
                   session: Session) -> None:
    if keyspace_exists(name,
                       session=session):
        logger.info(f'Dropping "{name}" keyspace.')
        drop_keyspace(name,
                      session=session)


if __name__ == '__main__':
    main()

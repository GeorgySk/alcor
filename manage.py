#!/usr/bin/env python3
import logging
import os
from contextlib import ExitStack

import click
from cassandra.cluster import (Cluster,
                               Session)
from cassandra.cqlengine import connection
from cassandra_helpers.connectable import \
    check_connection as check_cassandra_connection
from cassandra_helpers.keyspace import (keyspace_exists,
                                        create_keyspace,
                                        drop_keyspace)
from cassandra_helpers.models import sync_tables
from sqlalchemy.engine.url import make_url
from sqlalchemy.orm import sessionmaker
from sqlalchemy_helpers.connectable import (check_connection,
                                            get_engine,
                                            db_uri_to_str)
from sqlalchemy_utils import (database_exists,
                              create_database,
                              drop_database)

from alcor.cassandra_models import (CGroup,
                                    CStar,
                                    eliminations,
                                    luminosity_function,
                                    simulation,
                                    velocities,
                                    velocities_vs_magnitudes)
from alcor.config import PROJECT_NAME
from alcor.models.base import Base
from alcor.services.processing import run_processing
from alcor.services.simulations import run_simulations
from alcor.utils import load_settings

logger = logging.getLogger(__name__)


@click.group()
@click.pass_context
def main(ctx: click.Context) -> None:
    logging.basicConfig(level=logging.DEBUG)
    contact_points = [os.environ['CASSANDRA_RPC_ADDRESS']]
    port = int(os.environ['CASSANDRA_RPC_PORT'])
    ctx.obj = dict(contact_points=contact_points,
                   port=port,
                   db_uri=make_url(os.environ['POSTGRES_URI']))


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
              help='Cleans destination database.')
@click.pass_context
def simulate(ctx: click.Context,
             settings_path: str,
             project_dir: str,
             clean: bool) -> None:
    cluster_settings = ctx.obj

    db_uri = cluster_settings['db_uri']
    check_connection(db_uri)

    contact_points = cluster_settings['contact_points']
    port = cluster_settings['port']

    keyspace_name = PROJECT_NAME

    check_cassandra_connection(contact_points=contact_points,
                               port=port)
    with ExitStack() as stack:
        cluster = stack.enter_context(Cluster(contact_points=contact_points,
                                              port=port))
        engine = stack.enter_context(get_engine(db_uri))

        cassandra_session = cluster.connect()
        session_factory = sessionmaker(bind=engine)
        session = session_factory()

        if clean:
            clean_keyspace(keyspace_name,
                           session=cassandra_session)
            ctx.invoke(clean_db)

        init_cassandra(keyspace_name=keyspace_name,
                       session=cassandra_session)
        ctx.invoke(init_db)

        settings = load_settings(settings_path)
        os.chdir(project_dir)
        run_simulations(settings=settings,
                        cassandra_session=cassandra_session,
                        session=session)


@main.command()
@click.option('--filtration-method', '-m',
              type=click.Choice(['raw',
                                 'full',
                                 'restricted']),
              default='restricted',
              help='Raw data filtration method: '
                   '"raw" - do nothing, '
                   '"full" - only declination '
                   'and parallax selection criteria, '
                   '"restricted" - apply all criteria (default)')
@click.option('--nullify-radial-velocity', '-nrf',
              is_flag=True,
              help='Sets radial velocities to zero.')
@click.option('--w-luminosity-function', '-lf',
              is_flag=True,
              help='Prepare data for plotting luminosity function.')
@click.option('--w-velocities-clouds', '-uvw',
              is_flag=True,
              help='Prepare data for plotting velocity clouds.')
@click.option('--w-velocities-vs-magnitude', '-vm',
              is_flag=True,
              help='Prepare data for plots of velocities vs bol. magnitude.')
@click.option('--w-lepine-criterion', '-lcr',
              is_flag=True,
              help='Apply Lepine\'s criterion.')
@click.pass_context
def process(ctx: click.Context,
            filtration_method: str,
            nullify_radial_velocity: bool,
            w_luminosity_function: bool,
            w_velocities_clouds: bool,
            w_velocities_vs_magnitude: bool,
            w_lepine_criterion: bool) -> None:
    cluster_settings = ctx.obj

    db_uri = cluster_settings['db_uri']
    check_connection(db_uri)

    contact_points = cluster_settings['contact_points']
    port = cluster_settings['port']

    keyspace_name = PROJECT_NAME

    check_cassandra_connection(contact_points=contact_points,
                               port=port)
    with ExitStack() as stack:
        cluster = stack.enter_context(Cluster(contact_points=contact_points,
                                              port=port))
        engine = stack.enter_context(get_engine(db_uri))

        cassandra_session = cluster.connect()
        session_factory = sessionmaker(bind=engine)
        session = session_factory()

        init_cassandra(keyspace_name=keyspace_name,
                       session=cassandra_session)

        run_processing(filtration_method=filtration_method,
                       nullify_radial_velocity=nullify_radial_velocity,
                       w_luminosity_function=w_luminosity_function,
                       w_velocities_clouds=w_velocities_clouds,
                       w_velocities_vs_magnitude=w_velocities_vs_magnitude,
                       w_lepine_criterion=w_lepine_criterion,
                       session=session,
                       cassandra_session=cassandra_session)


@main.command(name='clean_db')
@click.pass_context
def clean_db(ctx: click.Context) -> None:
    """Removes Postgres database."""
    db_uri = ctx.obj['db_uri']
    db_uri_str = db_uri_to_str(db_uri)
    if database_exists(db_uri):
        logging.info(f'Cleaning "{db_uri_str}" database.')
        drop_database(db_uri)


@main.command(name='init_db')
@click.pass_context
def init_db(ctx: click.Context) -> None:
    """Creates Postgres database."""
    db_uri = ctx.obj['db_uri']
    db_uri_str = db_uri_to_str(db_uri)

    if not database_exists(db_uri):
        logging.info(f'Creating "{db_uri_str}" database.')
        create_database(db_uri)

    with get_engine(db_uri) as engine:
        logging.info(f'Creating "{db_uri_str}" database schema.')
        Base.metadata.create_all(bind=engine)


def init_cassandra(*,
                   keyspace_name: str,
                   session: Session) -> None:
    init_keyspace(keyspace_name,
                  session=session)
    session.set_keyspace(keyspace_name)
    connection.set_session(session)
    sync_tables(CGroup,
                CStar,
                simulation.CParameter,
                eliminations.CStarsCounter,
                luminosity_function.CPoint,
                velocities.CCloud,
                velocities.CLepineCaseUVCloud,
                velocities.CLepineCaseUWCloud,
                velocities.CLepineCaseVWCloud,
                velocities_vs_magnitudes.CCloud,
                velocities_vs_magnitudes.CLepineCaseUCloud,
                velocities_vs_magnitudes.CLepineCaseVCloud,
                velocities_vs_magnitudes.CLepineCaseWCloud,
                velocities_vs_magnitudes.CBin,
                velocities_vs_magnitudes.CLepineCaseUBin,
                velocities_vs_magnitudes.CLepineCaseVBin,
                velocities_vs_magnitudes.CLepineCaseWBin)


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

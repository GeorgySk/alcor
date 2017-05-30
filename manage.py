#!/usr/bin/env python3
import logging
import os

import click
from cassandra.cluster import (Cluster,
                               Session)
from cassandra.cqlengine import connection
from cassandra_helpers.connectable import check_connection
from cassandra_helpers.keyspace import (keyspace_exists,
                                        create_keyspace,
                                        drop_keyspace)
from cassandra_helpers.models import sync_tables

from alcor.config import PROJECT_NAME
from alcor.models import (Group,
                          Star,
                          eliminations,
                          luminosity_function,
                          simulation,
                          velocities,
                          velocities_vs_magnitudes)
from alcor.services.processing import run_processing
from alcor.services.plotting import draw_plots
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
def simulate(ctx: click.Context,
             settings_path: str,
             project_dir: str,
             clean: bool) -> None:
    cluster_settings = ctx.obj
    contact_points = cluster_settings['contact_points']
    port = cluster_settings['port']

    keyspace_name = PROJECT_NAME

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


@main.command()
@click.option('--filtration-method', '-m',
              type=click.Choice(['raw',
                                 'full',
                                 'restricted']),
              default='restricted',
              help='Raw data filtration method: '
                   '"raw" - do nothing, '
                   '"full" - only declination and parallax selection criteria, '
                   '"restricted" - apply all criteria (default)')
@click.option('--nullify-radial-velocity', '-nrf',
              is_flag=True,
              help='Sets radial velocities to zero.')
@click.option('--luminosity-function', '-lf',
              is_flag=True,
              help='Prepare data for plotting luminosity function.')
@click.option('--velocities-clouds', '-uvw',
              is_flag=True,
              help='Prepare data for plotting velocity clouds.')
@click.option('--velocities-vs-magnitude', '-vm',
              is_flag=True,
              help='Prepare data for plots of velocities vs bol. magnitude.')
@click.option('--lepine-criterion', '-lcr',
              is_flag=True,
              help='Apply Lepine\'s criterion.')
@click.pass_context
def process(ctx: click.Context,
            filtration_method: str,
            nullify_radial_velocity: bool,
            luminosity_function: bool,
            velocities_clouds: bool,
            velocities_vs_magnitude: bool,
            lepine_criterion: bool) -> None:
    cluster_settings = ctx.obj
    contact_points = cluster_settings['contact_points']
    port = cluster_settings['port']

    keyspace_name = PROJECT_NAME

    check_connection(contact_points=contact_points,
                     port=port)
    with Cluster(contact_points=contact_points,
                 port=port) as cluster:
        session = cluster.connect()
        init_db(keyspace_name=keyspace_name,
                session=session)

        run_processing(filtration_method=filtration_method,
                       nullify_radial_velocity=nullify_radial_velocity,
                       luminosity_function=luminosity_function,
                       velocities_clouds=velocities_clouds,
                       velocities_vs_magnitude=velocities_vs_magnitude,
                       lepine_criterion=lepine_criterion,
                       session=session)


@main.command()
@click.option('--luminosity-function', '-lf',
              is_flag=True,
              help='Plot luminosity function.')
@click.option('--velocities-vs-magnitude', '-vm',
              is_flag=True,
              help='Plot velocities vs bol. magnitude .')
@click.option('--velocity-clouds', '-uvw',
              is_flag=True,
              help='Plot velocity clouds.')
@click.option('--lepine-criterion', '-lcr',
              is_flag=True,
              help='Use data with applied Lepine\'s criterion.')
@click.pass_context
def plot(ctx: click.Context,
         luminosity_function: bool,
         velocities_vs_magnitude: bool,
         velocity_clouds: bool,
         lepine_criterion: bool) -> None:
    cluster_settings = ctx.obj
    contact_points = cluster_settings['contact_points']
    port = cluster_settings['port']
    keyspace_name = PROJECT_NAME
    check_connection(contact_points=contact_points,
                     port=port)
    with Cluster(contact_points=contact_points,
                 port=port) as cluster:
        session = cluster.connect()
        init_db(keyspace_name=keyspace_name,
                session=session)
        draw_plots(luminosity_function,
                   velocities_vs_magnitude,
                   velocity_clouds,
                   lepine_criterion)


def init_db(*,
            keyspace_name: str,
            session: Session) -> None:
    init_keyspace(keyspace_name,
                  session=session)
    session.set_keyspace(keyspace_name)
    connection.set_session(session)
    sync_tables(Group,
                Star,
                simulation.Parameter,
                eliminations.StarsCounter,
                luminosity_function.Point,
                velocities.Cloud,
                velocities.LepineCaseUVCloud,
                velocities.LepineCaseUWCloud,
                velocities.LepineCaseVWCloud,
                velocities_vs_magnitudes.Cloud,
                velocities_vs_magnitudes.LepineCaseUCloud,
                velocities_vs_magnitudes.LepineCaseVCloud,
                velocities_vs_magnitudes.LepineCaseWCloud,
                velocities_vs_magnitudes.Bin,
                velocities_vs_magnitudes.LepineCaseUBin,
                velocities_vs_magnitudes.LepineCaseVBin,
                velocities_vs_magnitudes.LepineCaseWBin)


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

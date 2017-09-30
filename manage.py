#!/usr/bin/env python3
import logging
import os
import uuid
from typing import Optional

import click
from sqlalchemy.engine.url import make_url
from sqlalchemy.orm import sessionmaker
from sqlalchemy_helpers.connectable import (check_connection,
                                            create_engine,
                                            db_uri_to_str)
from sqlalchemy_utils import (database_exists,
                              create_database,
                              drop_database)

from alcor.models.base import Base
from alcor.services import (simulations,
                            plots)
from alcor.services.common import FILTRATION_METHODS
from alcor.utils import load_settings

logger = logging.getLogger(__name__)


@click.group()
@click.pass_context
def main(ctx: click.Context) -> None:
    logging.basicConfig(format='%(filename)s %(funcName)s '
                               '%(levelname)s: %(message)s',
                        level=logging.DEBUG)
    ctx.obj = make_url(os.environ['POSTGRES_URI'])


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
    db_uri = ctx.obj
    check_connection(db_uri)

    with create_engine(db_uri) as engine:
        session_factory = sessionmaker(bind=engine)
        session = session_factory()

        if clean:
            ctx.invoke(clean_db)

        ctx.invoke(init_db)

        settings = load_settings(settings_path)
        os.chdir(project_dir)

        geometry = settings['geometry']
        precision = settings['precision']
        grid_parameters_settings = settings['grid']
        csv_parameters_settings = settings['csv']
        csv_parameters_info = {**csv_parameters_settings.get('common', {}),
                               **csv_parameters_settings.get(geometry, {})}
        grid_parameters_info = {**grid_parameters_settings.get('common', {}),
                                **grid_parameters_settings.get(geometry, {})}

        simulations.run(geometry=geometry,
                        precision=precision,
                        grid_parameters_info=grid_parameters_info,
                        csv_parameters_info=csv_parameters_info,
                        session=session)


@main.command()
# TODO: add other options
@click.option('--group_id',
              default=None,
              type=uuid.UUID,
              help='Process a group by id')
@click.option('--filtration-method', '-m',
              type=click.Choice(FILTRATION_METHODS),
              default='raw',
              help='Raw data filtration method: '
                   '"raw" - do nothing, '
                   '"full" - only declination '
                   'and parallax selection criteria, '
                   '"restricted" - apply all criteria (default)')
@click.option('--nullify-radial-velocity', '-nrv',
              is_flag=True,
              help='Sets radial velocities to zero.')
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
@click.option('--heatmap', '-hm',
              type=click.Choice(['velocities',
                                 'coordinates']),
              help='Plot heatmaps for: '
                   '"velocities" - velocities, '
                   '"coordinates" - coordinates')
@click.option('--toomre-diagram', '-toomre',
              is_flag=True,
              help='Plot Toomre diagram.')
@click.option('--ugriz-color-color-diagram', '-ugriz',
              is_flag=True,
              help='Plot color-color diagrams for ugriz photometry.')
@click.pass_context
def plot(ctx: click.Context,
         group_id: Optional[uuid.UUID],
         filtration_method: str,
         nullify_radial_velocity: bool,
         luminosity_function: bool,
         velocities_vs_magnitude: bool,
         velocity_clouds: bool,
         lepine_criterion: bool,
         heatmap: str,
         toomre_diagram: bool,
         ugriz_color_color_diagram: bool) -> None:
    db_uri = ctx.obj
    check_connection(db_uri)

    with create_engine(db_uri) as engine:
        session_factory = sessionmaker(bind=engine)
        session = session_factory()
        plots.draw(group_id,
                   filtration_method,
                   nullify_radial_velocity,
                   luminosity_function,
                   velocities_vs_magnitude,
                   velocity_clouds,
                   lepine_criterion,
                   heatmap,
                   toomre_diagram,
                   ugriz_color_color_diagram,
                   session=session)


@main.command(name='clean_db')
@click.pass_context
def clean_db(ctx: click.Context) -> None:
    """Removes Postgres database."""
    db_uri = ctx.obj
    db_uri_str = db_uri_to_str(db_uri)
    if database_exists(db_uri):
        logging.info(f'Cleaning "{db_uri_str}" database.')
        drop_database(db_uri)


@main.command(name='init_db')
@click.pass_context
def init_db(ctx: click.Context) -> None:
    """Creates Postgres database."""
    db_uri = ctx.obj
    db_uri_str = db_uri_to_str(db_uri)

    if not database_exists(db_uri):
        logging.info(f'Creating "{db_uri_str}" database.')
        create_database(db_uri)

    with create_engine(db_uri) as engine:
        logging.info(f'Creating "{db_uri_str}" database schema.')
        Base.metadata.create_all(bind=engine)


if __name__ == '__main__':
    main()

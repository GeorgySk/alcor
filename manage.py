#!/usr/bin/env python3
import logging
import os

import click
from alcor.services import simulations
from alcor.utils import load_settings

logger = logging.getLogger(__name__)


@click.group()
def main() -> None:
    logging.basicConfig(format='%(filename)s %(funcName)s '
                               '%(levelname)s: %(message)s',
                        level=logging.DEBUG)


@main.command()
@click.option('--settings-path', '-p',
              default='settings.yml',
              type=click.Path(),
              help='Settings file path '
                   '(absolute or relative, '
                   'default "settings.yml").')
@click.option('--project-dir', '-d',
              required=True,
              default='test_project',
              type=click.Path(),
              help='Project directory path '
                   '(absolute or relative).')
def simulate(settings_path: str,
             project_dir: str) -> None:
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
                    csv_parameters_info=csv_parameters_info)


if __name__ == '__main__':
    main()

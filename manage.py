#!/usr/bin/env python3
import os
import uuid
from itertools import product
from subprocess import call
from typing import (Union,
                    Iterable,
                    Dict,
                    Tuple,
                    List)

import click
import logging

from alcor.utils import load_settings

logger = logging.getLogger(__name__)

NumericType = Union[int, float]


@click.group()
def main() -> None:
    logging.basicConfig(level=logging.DEBUG)


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
def run(settings_path: str,
        project_dir: str) -> None:
    settings = load_settings(settings_path)
    model_type = settings['model_type']
    precision = settings['precision']
    parameters_info = settings['parameters']
    os.chdir(project_dir)
    max_output_file_name_length = 5
    output_file_extension = '.res'
    for parameters_values in generate_parameters_values(
            parameters_info=parameters_info,
            precision=precision):
        output_file_name = (uuid.uuid4().hex[:max_output_file_name_length]
                            + output_file_extension)
        simulate(parameters_values=parameters_values,
                 model_type=model_type,
                 output_file_name=output_file_name)


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
    call(args)


def generate_parameters_values(
        parameters_info: Dict[str, Dict[str, NumericType]],
        precision: int) -> Iterable[Dict[str, NumericType]]:
    parameters_values_ranges_by_names = dict(
        generate_parameters_values_ranges_by_names(
            parameters_info=parameters_info,
            precision=precision))
    parameters_names = list(parameters_values_ranges_by_names.keys())
    for values in product(*parameters_values_ranges_by_names.values()):
        yield dict(zip(parameters_names,
                       values))


def generate_parameters_values_ranges_by_names(
        *,
        parameters_info: Dict[str, Dict[str, NumericType]],
        precision: int) -> Iterable[Tuple[str, List[NumericType]]]:
    for parameter_name, parameter_settings in parameters_info.items():
        start_value = round(parameter_settings['start'], precision)
        step_size = round(parameter_settings['step'], precision)
        values_count = parameter_settings['count']
        values_range = [
            round(start_value + value_number * step_size,
                  precision)
            for value_number in range(values_count)
            ]
        yield parameter_name, values_range


if __name__ == '__main__':
    main()

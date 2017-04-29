#!/usr/bin/env python3
from itertools import product
from typing import (Union,
                    Iterable,
                    Dict,
                    Tuple,
                    List)

import click

from alcor.utils import load_settings

NumericType = Union[int, float]


@click.group()
def main() -> None:
    pass


@main.command()
@click.option('--settings-path', '-p',
              default='settings.yml',
              type=click.Path(),
              help='Settings file path '
                   '(absolute or relative, '
                   'default "settings.yml").')
def run(settings_path: str) -> None:
    settings = load_settings(settings_path)
    precision = settings['precision']
    parameters_info = settings['parameters']
    parameters_values = list(generate_parameters_values(
        parameters_info=parameters_info,
        precision=precision))


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

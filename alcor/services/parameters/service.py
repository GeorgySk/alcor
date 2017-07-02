import csv
import logging
from itertools import product
from typing import (Iterable,
                    Dict,
                    Tuple,
                    List)

from alcor.types import NumericType


logger = logging.getLogger(__name__)


def generate_parameters_values(*,
                               parameters: Dict[str, Dict[str, NumericType]],
                               precision: int,
                               geometry: str
                               ) -> Iterable[Dict[str, NumericType]]:
    parameters_values_ranges_by_names = dict(
        generate_parameters_values_ranges_by_names(parameters_dict=parameters,
                                                   precision=precision,
                                                   geometry=geometry))

    parameters_keys_read_from_file = []

    for key, value in parameters_values_ranges_by_names.items():
        try:
            if isinstance(value[0], str):
                parameters_keys_read_from_file.append(key)
        except IndexError:
            print(f'No values found in the file for parameter {key}')

    if parameters_keys_read_from_file:
        file_lines_count = len(parameters_values_ranges_by_names[
                                   parameters_keys_read_from_file[0]])

        all_simulations_parameters = []

        for file_line in range(file_lines_count):
            one_simulation_parameters = [
                parameters_values_ranges_by_names[key][file_line]
                if key in parameters_keys_read_from_file
                else parameters_values_ranges_by_names[key]
                for key, value in parameters_values_ranges_by_names.items()]
            all_simulations_parameters.append(one_simulation_parameters)

        product_list = [tuple(product(*one_simulation_parameters))
                        for one_simulation_parameters in
                        all_simulations_parameters]

        unpacked_product_list = [tuple_inner_instance
                                 for tuple_instance in product_list
                                 for tuple_inner_instance in tuple_instance]

        parameters_names = list(parameters_values_ranges_by_names.keys())

        for values in unpacked_product_list:
            yield dict(zip(parameters_names, values))
    else:
        parameters_names = list(parameters_values_ranges_by_names.keys())

        for values in product(*parameters_values_ranges_by_names.values()):
            yield dict(zip(parameters_names,
                           values))


def generate_parameters_values_ranges_by_names(
        *,
        parameters_dict: Dict[str, Dict[str, NumericType]],
        precision: int,
        geometry: str) -> Iterable[Tuple[str, List[NumericType]]]:

    variable_simulations_parameters = parameters_dict['commons']

    if geometry == 'cone':
        variable_simulations_parameters.update(parameters_dict['cone'])

    for parameter_name, parameter_settings in (variable_simulations_parameters.
                                               items()):

        # Parameter is const during all simulations:
        if isinstance(parameter_settings, float):
            values_range = [round(parameter_settings, precision)]
            yield parameter_name, values_range

        elif isinstance(parameter_settings, dict):

            # Parameter is generated for every simulation:
            if all(key in parameter_settings
                   for key in ['start', 'step', 'count']):
                start_value = round(parameter_settings['start'], precision)
                step_size = round(parameter_settings['step'], precision)
                values_count = parameter_settings['count']
                values_range = [
                    round(start_value + value_number * step_size, precision)
                    for value_number in range(values_count)]
                yield parameter_name, values_range

            # Parameter is read from csv file:
            elif all(key in parameter_settings
                     for key in ['csv', 'column']):
                file_path = parameter_settings['csv']
                column = parameter_settings['column']
                values_range = []
                with open(file_path, 'r') as file:
                    file_reader = csv.reader(file,
                                             delimiter=' ',
                                             skipinitialspace=True)
                    for row in file_reader:
                        try:
                            values_range.append(row[column - 1])
                        except IndexError:
                            logger.error('Inconsistent number of lines in '
                                         'csv file with parameters')
                yield parameter_name, values_range

            else:
                logger.error('Wrong set of keys. Possible sets are '
                             '(start, step, count) and (csv, column)')
        else:
            logger.error(f'Wrong type of settings {parameter_settings}. '
                         f'Expected float or dict, got '
                         f'{type(parameter_settings)} instead')

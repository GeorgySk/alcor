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

    # Parameters read from file are saved in lists of strings
    str_lists_indexes = []

    for list_index, list_instance in enumerate(
            parameters_values_ranges_by_names.values()):
        if type(list_instance[0]) is str:
            str_lists_indexes.append(list_index)

    file_lines_count = len(list(
        parameters_values_ranges_by_names.values())[str_lists_indexes[0]])

    all_simulations_parameters = []
    one_simulation_parameters = []
    for line_index in range(file_lines_count):
        for list_index, list_instance in enumerate(
                parameters_values_ranges_by_names.values()):
            if list_index in str_lists_indexes:
                one_simulation_parameters.append([list_instance[line_index]])
            else:
                one_simulation_parameters.append(list_instance)
        all_simulations_parameters.append(one_simulation_parameters)
        one_simulation_parameters = []

    product_list = []
    for one_simulation_parameters in all_simulations_parameters:
        product_list.append(tuple(product(*one_simulation_parameters)))

    unpacked_product_list = []
    for tuple_instance in product_list:
        for tuple_inner_instance in tuple_instance:
            unpacked_product_list.append(tuple_inner_instance)

    parameters_names = list(parameters_values_ranges_by_names.keys())

    for values in unpacked_product_list:
        yield dict(zip(parameters_names, values))


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
        if type(parameter_settings) is float:
            values_range = [round(parameter_settings, precision)]
            yield parameter_name, values_range

        elif type(parameter_settings) is dict:

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
                        values_range.append(row[column - 1])
                yield parameter_name, values_range

            else:
                logger.error('Wrong set of keys. Possible sets are '
                             '(start, step, count) and (csv, column)')
        else:
            logger.error(f'Wrong type of settings {parameter_settings}. '
                         f'Expected float or dict, got '
                         f'{type(parameter_settings)} instead')

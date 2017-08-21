import csv
import logging
from itertools import product
from typing import (Iterable,
                    Dict,
                    Tuple,
                    List)

from alcor.types import NumericType

GRID_SETTINGS_KEYS = ['start', 'step', 'count']

logger = logging.getLogger(__name__)


def generate_parameters_values(
        parameters_info: Dict[str, Dict[str, NumericType]],
        precision: int,
        geometry: str) -> Iterable[Dict[str, NumericType]]:
    current_geometry_parameters_info = {**parameters_info['commons'],
                                        **parameters_info[geometry]}

    variable_parameters_info = {}
    non_variable_parameters_info = {}
    for parameter, value in (current_geometry_parameters_info.items()):
        is_variable_parameter = (isinstance(value, dict)
                                 and all(key in value
                                         for key in GRID_SETTINGS_KEYS))
        if is_variable_parameter:
            variable_parameters_info[parameter] = value
        else:
            non_variable_parameters_info[parameter] = value

    parameters_values_ranges_by_names = dict(
        generate_parameters_values_ranges_by_names(
            parameters_info=variable_parameters_info,
            precision=precision))
    parameters_names = list(parameters_values_ranges_by_names.keys())
    for values in product(*parameters_values_ranges_by_names.values()):
        variable_parameters_dict = dict(zip(parameters_names, values))
        output_dict = {**variable_parameters_dict,
                       **non_variable_parameters_info}
        yield output_dict


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

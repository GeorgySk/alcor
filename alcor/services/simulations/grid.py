import logging
from itertools import product
from typing import (Iterator,
                    Dict,
                    Tuple,
                    List)

from alcor.types import (NumericType,
                         GridParametersInfoType)

logger = logging.getLogger(__name__)


def parameters_values(
        parameters_info: GridParametersInfoType,
        precision: int) -> Iterator[Dict[str, NumericType]]:
    values_ranges_by_names = dict(
            parameters_names_values_ranges(
                    parameters_info=parameters_info,
                    precision=precision))
    parameters_names = list(values_ranges_by_names)
    for values in product(*values_ranges_by_names.values()):
        yield dict(zip(parameters_names, values))


def parameters_names_values_ranges(
        *,
        parameters_info: GridParametersInfoType,
        precision: int) -> Iterator[Tuple[str, List[NumericType]]]:
    for name, value_or_grid in parameters_info.items():
        try:
            start_value = round(value_or_grid['start'], precision)
        except TypeError:
            # single value
            yield name, [round(value_or_grid, precision)]
        else:
            # grid
            step_size = round(value_or_grid['step'], precision)
            values_count = value_or_grid['count']
            values_range = [
                round(start_value + value_number * step_size,
                      precision)
                for value_number in range(values_count)
            ]
            yield name, values_range

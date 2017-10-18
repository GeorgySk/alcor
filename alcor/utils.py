import logging
import operator
from collections import OrderedDict
from functools import reduce
from typing import (Any,
                    Union,
                    Hashable,
                    Iterable,
                    Iterator,
                    Container,
                    Mapping,
                    Dict,
                    Tuple)

import yaml

from alcor.models import (STAR_PARAMETERS_NAMES,
                          Group,
                          Star)

logger = logging.getLogger(__name__)


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


def join_str(items: Iterable[Any],
             sep: str = ', ') -> str:
    return sep.join(map(str, items))


def parse_stars(lines: Iterator[str],
                *,
                group: Group,
                possible_columns_names: Container[str] = STAR_PARAMETERS_NAMES
                ) -> Iterator[Star]:
    header = next(lines).split()
    validate_header(header,
                    possible_columns_names=possible_columns_names)

    for line in lines:
        parts = line.split()
        params = map(str_to_float, parts)
        values = OrderedDict(zip(header,
                                 params))
        yield Star(group_id=group.id,
                   **values)


def validate_header(header: Iterable[str],
                    *,
                    possible_columns_names: Container[str]) -> None:
    unknown_columns_names = [column_name
                             for column_name in header
                             if column_name not in possible_columns_names]
    if unknown_columns_names:
        err_msg = ('Unknown columns names: "{columns_names}".'
                   .format(columns_names='", "'.join(unknown_columns_names)))
        raise ValueError(err_msg)


def str_to_float(string: str) -> Union[str, float]:
    try:
        return float(string)
    except ValueError:
        return string


def zip_mappings(*mappings: Mapping[Hashable, Any]
                 ) -> Iterator[Tuple[Hashable, Tuple[Any, ...]]]:
    keys_sets = map(set, mappings)
    common_keys = reduce(set.intersection, keys_sets)
    for key in common_keys:
        yield key, tuple(map(operator.itemgetter(key), mappings))

import logging
import operator
from collections import OrderedDict
from functools import reduce
from typing import (Any,
                    Union,
                    Hashable,
                    Iterator,
                    Container,
                    Mapping,
                    Dict,
                    Tuple,
                    List)

import yaml

from alcor.models import (Group,
                          Star)

logger = logging.getLogger(__name__)


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


def parse_stars(lines: Iterator[str],
                *,
                group: Group,
                columns_names: List[str]) -> Iterator[Star]:
    group_id = group.id
    for line in lines:
        parts = line.split()
        params = map(str_to_float, parts)
        values = OrderedDict(zip(columns_names,
                                 params))
        yield Star(group_id=group_id,
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

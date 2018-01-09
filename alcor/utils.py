import logging
import operator
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

logger = logging.getLogger(__name__)


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


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

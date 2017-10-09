import logging
import operator
from collections import OrderedDict
from functools import reduce
from typing import (Any,
                    Union,
                    Iterable,
                    Iterator,
                    Mapping,
                    Hashable,
                    Dict)

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
                group: Group) -> Iterator[Star]:
    headers = next(lines).split()
    for header in headers:
        if not (header in STAR_PARAMETERS_NAMES):
            logger.error(f'There is no parameter {header} in '
                         f'STAR_PARAMETERS_NAMES')
    for line in lines:
        parts = line.split()
        params = map(str_to_float, parts)
        values = OrderedDict(zip(headers,
                                 params))
        yield Star(group_id=group.id,
                   **values)


def str_to_float(string: str) -> Union[str, float]:
    try:
        return float(string)
    except ValueError:
        return string


def zip_mappings(*mappings: Mapping[Hashable, Any]):
    keys_sets = map(set, mappings)
    common_keys = reduce(set.intersection, keys_sets)
    for key in common_keys:
        yield key, tuple(map(operator.itemgetter(key), mappings))

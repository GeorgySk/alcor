import decimal
import logging
from collections import OrderedDict
from decimal import Decimal
from typing import (Any,
                    Iterable,
                    Iterator,
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
        params = map(str_to_decimal, parts)
        values = OrderedDict(zip(headers,
                                 params))
        yield Star(group_id=group.id,
                   **values)


def str_to_decimal(string: str) -> Any:
    try:
        return Decimal(string)
    except decimal.InvalidOperation:
        return string

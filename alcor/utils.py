import logging
from collections import OrderedDict
from decimal import Decimal
from typing import (Any,
                    Iterable,
                    Dict,
                    Tuple)

import yaml

from alcor.models import (Group,
                          Star)
from alcor.models import STAR_PARAMETERS_NAMES
from alcor.types import ColumnValueType

logger = logging.getLogger(__name__)


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


def join_str(items: Iterable[Any],
             sep: str = ', ') -> str:
    return sep.join(map(str, items))


def parse_stars(lines: Iterable[str],
                group: Group
                ) -> Iterable[Star]:
    headers = next(lines).split()
    for header in headers:
        if not (header in STAR_PARAMETERS_NAMES):
            logger.error(f'There is no parameter {header} in '
                         f'STAR_PARAMETERS_NAMES')
    for line in lines:
        parts = line.split()
        params = map(Decimal, parts)
        values = OrderedDict(zip(headers,
                                 params))
        yield Star(group_id=group.id,
                   **values)


def get_columns(rows: Iterable[Tuple[str, ...]],
                data_type: ColumnValueType = float
                ) -> Iterable[Tuple[ColumnValueType, ...]]:
    converted_rows = (map(data_type, row)
                      for row in rows)
    return zip(*converted_rows)

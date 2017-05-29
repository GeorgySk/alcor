import uuid
from collections import OrderedDict
from typing import (Any,
                    Iterable,
                    Dict)

import yaml
from decimal import Decimal

from alcor.models import Star
from alcor.models.star import STAR_PARAMETERS_NAMES


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


def join_str(items: Iterable[Any], sep: str = ', ') -> str:
    return sep.join(map(str, items))


def parse_stars(lines: Iterable[str],
                group_id: uuid.UUID
                ) -> Iterable[Star]:
    for line in lines:
        parts = line.split()
        params = map(Decimal, parts)
        values = OrderedDict(zip(STAR_PARAMETERS_NAMES,
                                 params))
        yield Star(group_id=group_id,
                   **values)

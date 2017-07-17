from collections import OrderedDict
from decimal import Decimal
from itertools import chain
from typing import (Any,
                    Iterable,
                    Dict,
                    Tuple)

import yaml

from alcor.cassandra_models import (STAR_PARAMETERS_NAMES,
                                    CGroup,
                                    CStar)
from alcor.models import (Group,
                          Star)


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


def join_str(items: Iterable[Any], sep: str = ', ') -> str:
    return sep.join(map(str, items))


def parse_stars(lines: Iterable[str],
                c_group: CGroup,
                group: Group
                ) -> Iterable[Tuple[CStar, Star]]:
    for line in lines:
        parts = line.split()
        params = chain(map(Decimal, parts[:-1]),
                       # spectral type is integer
                       [int(parts[-1])])
        values = OrderedDict(zip(STAR_PARAMETERS_NAMES,
                                 params))
        yield (CStar(group_id=c_group.id,
                     **values),
               Star(group_id=group.id,
                    **values))

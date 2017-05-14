from typing import (Any,
                    Iterable,
                    Dict)

import yaml


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)


def join_str(items: Iterable[Any], sep: str = ', ') -> str:
    return sep.join(map(str, items))

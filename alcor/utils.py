from typing import (Any,
                    Dict)

import yaml


def load_settings(path: str
                  ) -> Dict[str, Any]:
    with open(path) as file:
        return yaml.safe_load(file)

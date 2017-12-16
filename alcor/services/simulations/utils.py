from typing import (Any,
                    Iterable)

import numpy as np


def immutable_array(elements: Iterable[Any]) -> np.ndarray:
    result = np.array(elements)
    result.setflags(write=False)
    return result


def linear_function(values: np.ndarray,
                    *,
                    factor: float,
                    const: float) -> np.ndarray:
    return values * factor + const

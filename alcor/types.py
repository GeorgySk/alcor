from typing import (Union,
                    Dict,
                    Callable)

import numpy as np

NumericType = Union[int, float]
GridParametersInfoType = Dict[str, Union[NumericType,
                                         Dict[str, NumericType]]]
CSVParametersInfoType = Dict[str, Dict[str, Union[str, int]]]
ArrayOperatorType = Callable[[np.ndarray], np.ndarray]

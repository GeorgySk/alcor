from typing import (Union,
                    Callable,
                    Iterable,
                    Dict,
                    Tuple)

import numpy as np

NumericType = Union[int, float]
GridParametersInfoType = Dict[str, Union[NumericType,
                                         Dict[str, NumericType]]]
CSVParametersInfoType = Dict[str, Dict[str, Union[str, int]]]
ParametersValuesType = Dict[str, Union[NumericType,
                                       Dict[str, Union[str, NumericType]]]]
GaussianGeneratorType = Callable[[Union[float, np.ndarray, Iterable],
                                  Union[float, np.ndarray, Iterable],
                                  Union[int, Iterable[int], None]],
                                 Union[np.ndarray, int, float, complex]]
GeneratorType = Callable[[float, float, float], np.ndarray]
UnitRangeGeneratorType = Callable[[Tuple[int, ...]], np.ndarray]
ArrayOperatorType = Callable[[np.ndarray], np.ndarray]

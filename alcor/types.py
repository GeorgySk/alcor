from typing import (Union,
                    Dict)

NumericType = Union[int, float]
GridParametersInfoType = Dict[str, Union[NumericType,
                                         Dict[str, NumericType]]]
CSVParametersInfoType = Dict[str, Dict[str, Union[str, int]]]
ParametersValuesType = Dict[str, Union[NumericType,
                                       Dict[str, Union[str, NumericType]]]]
CoolingSequenceType = Dict[float, Dict[str, np.ndarray]]

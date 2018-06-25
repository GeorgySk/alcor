from typing import (Union,
                    Dict)

NumericType = Union[int, float]
GridParametersInfoType = Dict[str, Union[NumericType,
                                         Dict[str, NumericType]]]
CSVParametersInfoType = Dict[str, Dict[str, Union[str, int]]]

import uuid
from datetime import datetime
from numbers import Number
from typing import (Union,
                    Callable,
                    Dict,
                    Tuple,
                    List)

from alcor.models import Star

NumericType = Union[int, float]
StarsBinsType = List[List[Star]]
RowType = Tuple[Number, ...]
ColumnValueType = Union[int,
                        bool,
                        float,
                        str,
                        datetime,
                        uuid.UUID,
                        None]
BolometricIndexType = Callable[[float], int]
StarBolometricIndexType = Callable[[Star], int]
GridParametersInfoType = Dict[str, Union[NumericType,
                                         Dict[str, NumericType]]]
CSVParametersInfoType = Dict[str, Dict[str, Union[str, int]]]

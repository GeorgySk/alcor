import uuid
from datetime import datetime
from numbers import Number
from typing import (Union,
                    Callable,
                    Dict,
                    Tuple,
                    List)

import numpy as np

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
ParametersValuesType = Dict[str, Union[NumericType,
                                       Dict[str, Union[str, NumericType]]]]
BolometricIndexType = Callable[[float], int]
StarBolometricIndexType = Callable[[Star], int]
CoolingSequencesType = Dict[int, Dict[str, np.ndarray]]

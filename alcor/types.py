import uuid
from datetime import datetime
from numbers import Number
from typing import (Union,
                    List,
                    Tuple)

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

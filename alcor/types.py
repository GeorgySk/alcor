from datetime import datetime
from numbers import Number
from typing import (Union,
                    Optional,
                    Callable,
                    Dict,
                    Tuple,
                    List, Iterable)

from cassandra.query import Statement

from alcor.cassandra_models import CStar

NumericType = Union[int, float]
StarsBinsType = List[List[CStar]]
RowType = Tuple[Number, ...]
ColumnValueType = Union[int, bool,
                        float, str,
                        datetime, None]
RecordType = Dict[str, ColumnValueType]
# callback argument will be `None` for "insert" function
ResponseType = Optional[List[RecordType]]
CallbackType = Callable[[ResponseType], None]
StatementType = Union[str, Statement]
StatementParametersType = Iterable[ColumnValueType]

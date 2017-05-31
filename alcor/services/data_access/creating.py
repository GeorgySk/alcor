from typing import Iterable

from cassandra.cluster import Session
from cassandra.cqlengine.models import Model

from alcor.types import (CallbackType,
                         StatementType)
from .callbacks import (empty_callback,
                        add_callback)


def insert(*,
           instances: Iterable[Model],
           statement: StatementType,
           session: Session,
           callback: CallbackType = empty_callback
           ) -> None:
    prepared_statement = session.prepare(statement)
    for instance in instances:
        future = session.execute_async(prepared_statement,
                                       parameters=instance.values())
        add_callback(future=future,
                     callback=callback)

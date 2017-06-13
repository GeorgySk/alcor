from typing import Iterable

from cassandra.cluster import Session
from cassandra.cqlengine.models import Model

from alcor.types import CallbackType
from .callbacks import empty_callback
from .execution import execute_prepared_statement


def insert(*,
           statement: str,
           instances: Iterable[Model],
           session: Session,
           callback: CallbackType = empty_callback
           ) -> None:
    parameters_collection = [instance.values()
                             for instance in instances]
    execute_prepared_statement(statement=statement,
                               parameters_collection=parameters_collection,
                               callback=callback,
                               session=session)

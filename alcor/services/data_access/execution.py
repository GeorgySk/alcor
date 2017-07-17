from typing import (Optional,
                    Dict,
                    List, Iterable)

from cassandra.cluster import Session
from cassandra.cqlengine.statements import BaseCQLStatement
from cassandra.query import SimpleStatement

from alcor.types import (CallbackType,
                         StatementType,
                         RecordType,
                         ColumnValueType,
                         StatementParametersType)
from .callbacks import (add_callback,
                        empty_callback)


def execute_base_statement(*,
                           base_statement: BaseCQLStatement,
                           session: Session,
                           callback: CallbackType = None,
                           statement_cls=SimpleStatement,
                           **kwargs
                           ) -> Optional[List[RecordType]]:
    parameters = base_statement.get_context()
    statement = statement_cls(str(base_statement),
                              **kwargs)
    return execute_statement(statement=statement,
                             parameters=parameters,
                             session=session,
                             callback=callback)


def execute_statement(*,
                      statement: StatementType,
                      parameters: Dict[str, ColumnValueType] = None,
                      session: Session,
                      callback: CallbackType = None
                      ) -> Optional[List[RecordType]]:
    future = session.execute_async(statement,
                                   parameters=parameters)
    return add_callback(future=future,
                        callback=callback)


def execute_prepared_statement(*,
                               statement: str,
                               parameters_collection: Iterable[StatementParametersType],
                               callback: CallbackType = empty_callback,
                               session: Session) -> None:
    prepared_statement = session.prepare(statement)
    for parameters in parameters_collection:
        future = session.execute_async(prepared_statement,
                                       parameters=parameters)
        add_callback(future=future,
                     callback=callback)

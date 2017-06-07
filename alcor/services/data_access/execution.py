from typing import (Optional,
                    Dict,
                    List)

from cassandra.cluster import Session
from cassandra.cqlengine.statements import BaseCQLStatement
from cassandra.query import SimpleStatement

from alcor.services.data_access.callbacks import add_callback
from alcor.types import (CallbackType,
                         StatementType,
                         RecordType,
                         ColumnValueType)


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

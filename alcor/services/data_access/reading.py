from typing import (Optional,
                    Dict,
                    List)

from cassandra.cluster import Session
from cassandra.cqlengine.query import AbstractQuerySet
from cassandra.cqlengine.statements import SelectStatement
from cassandra.query import SimpleStatement

from alcor.types import (StatementType,
                         CallbackType,
                         RecordType,
                         ColumnValueType)
from .callbacks import add_callback
from .statements import query_to_select_statement


def fetch_query(*,
                query: AbstractQuerySet,
                session: Session,
                callback: CallbackType = None,
                statement_cls=SimpleStatement,
                **kwargs
                ) -> Optional[List[RecordType]]:
    select_statement = query_to_select_statement(query)
    return fetch_select_statement(select_statement=select_statement,
                                  session=session,
                                  callback=callback,
                                  statement_cls=statement_cls,
                                  **kwargs)


def fetch_select_statement(*,
                           select_statement: SelectStatement,
                           session: Session,
                           callback: CallbackType = None,
                           statement_cls=SimpleStatement,
                           **kwargs
                           ) -> Optional[List[RecordType]]:
    parameters = select_statement.get_context()
    statement = statement_cls(str(select_statement),
                              **kwargs)
    return fetch_statement(statement=statement,
                           parameters=parameters,
                           session=session,
                           callback=callback)


def fetch_statement(*,
                    statement: StatementType,
                    parameters: Dict[str, ColumnValueType] = None,
                    session: Session,
                    callback: CallbackType = None
                    ) -> Optional[List[RecordType]]:
    future = session.execute_async(statement,
                                   parameters=parameters)
    return add_callback(future=future,
                        callback=callback)

from typing import (Optional,
                    List)

from cassandra.cluster import Session
from cassandra.cqlengine.query import AbstractQuerySet
from cassandra.query import SimpleStatement

from alcor.types import (CallbackType,
                         RecordType)
from .execution import execute_base_statement
from .statements import query_to_select_statement


def fetch(*,
          query: AbstractQuerySet,
          session: Session,
          callback: CallbackType = None,
          statement_cls=SimpleStatement,
          **kwargs
          ) -> Optional[List[RecordType]]:
    select_statement = query_to_select_statement(query)
    return execute_base_statement(base_statement=select_statement,
                                  session=session,
                                  callback=callback,
                                  statement_cls=statement_cls,
                                  **kwargs)

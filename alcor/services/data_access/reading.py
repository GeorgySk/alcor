from typing import (Optional,
                    List)

from cassandra.cluster import Session

from alcor.types import (StatementType,
                         CallbackType,
                         RecordType)
from .callbacks import add_callback


def fetch(*,
          statement: StatementType,
          session: Session,
          callback: CallbackType = None
          ) -> Optional[List[RecordType]]:
    future = session.execute_async(statement)
    return add_callback(future=future,
                        callback=callback)

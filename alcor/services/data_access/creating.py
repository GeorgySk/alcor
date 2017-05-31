import logging
from typing import (Iterable,
                    List)

from cassandra.cluster import Session
from cassandra.cqlengine.models import Model

from alcor.types import CallbackType
from alcor.utils import join_str
from .callbacks import (empty_callback,
                        add_callback)

logger = logging.getLogger(__name__)
BIND_ALIAS = '?'


def insert(*,
           instances: Iterable[Model],
           session: Session,
           callback: CallbackType = empty_callback
           ) -> None:
    instances = iter(instances)

    try:
        instance = next(instances)
    except StopIteration:
        return
    table_name = instance.__table_name__

    query = insert_query(table_name=table_name,
                         columns_names=instance.keys())
    prepared_statement = session.prepare(query)
    future = session.execute_async(prepared_statement,
                                   parameters=instance.values())
    add_callback(future=future,
                 callback=callback)
    for num, instance in enumerate(instances, start=2):
        future = session.execute_async(prepared_statement,
                                       parameters=instance.values())
        add_callback(future=future,
                     callback=callback)
    msg = (f'Successfully finished processing '
           f'"{table_name}" table\'s records, '
           f'number of instances processed: {num}.')
    logger.debug(msg)


def insert_query(*,
                 table_name: str,
                 columns_names: List[str]) -> str:
    binds_aliases = [BIND_ALIAS] * len(columns_names)
    columns_names_str = join_str(columns_names)
    binds_aliases_str = join_str(binds_aliases)
    return (f'INSERT INTO {table_name} ({columns_names_str}) '
            f'VALUES ({binds_aliases_str})')

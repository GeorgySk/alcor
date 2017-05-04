import logging
from functools import partial
from itertools import islice
from typing import (Iterable,
                    List)

from cassandra.cluster import Session
from cassandra.cqlengine.models import Model
from cassandra.query import (BatchStatement,
                             SimpleStatement)

from alcor.utils import join_str

logger = logging.getLogger(__name__)

BIND_ALIAS = '%s'
PREPARED_BIND_ALIAS = '?'

MAX_STATEMENTS_IN_BATCH = 65_535


def insert_by_batches(
        *,
        instances: Iterable[Model],
        table_name: str,
        session: Session,
        batch_size: int) -> None:
    validate_batch_size(batch_size)

    instances = iter(instances)
    while True:
        batch_instances = list(islice(instances,
                                      batch_size))
        if not batch_instances:
            return
        insert_batch(instances=batch_instances,
                     table_name=table_name,
                     session=session)


def insert_batch(*,
                 instances: Iterable[Model],
                 table_name: str,
                 session: Session) -> None:
    query = partial(insert_query,
                    table_name)
    batch = BatchStatement()
    for num, instance in enumerate(instances, start=1):
        query_string = query(instance.keys())
        batch.add(SimpleStatement(query_string),
                  instance.values())

    def callback(resp: None) -> None:
        msg = (f'Successfully processed '
               f'"{table_name}" table\'s instances, '
               f'number of instances processed: {num}.')
        logger.debug(msg)

    future = session.execute_async(batch)
    future.add_callback(callback)


def insert_query(table_name: str,
                 columns_names: List[str]) -> str:
    binds_aliases = [BIND_ALIAS] * len(columns_names)
    columns_names_str = join_str(columns_names)
    binds_aliases_str = join_str(binds_aliases)
    return (f'INSERT INTO {table_name} ({columns_names_str}) '
            f'VALUES ({binds_aliases_str})')


def validate_batch_size(size: int) -> None:
    if size > MAX_STATEMENTS_IN_BATCH:
        err_msg = ('Batch statement cannot contain '
                   f'more than {MAX_STATEMENTS_IN_BATCH} statements.')
        raise ValueError(err_msg)


def insert(*,
           instances: Iterable[Model],
           session: Session) -> None:
    instances = iter(instances)

    try:
        instance = next(instances)
    except StopIteration:
        return
    table_name = instance.__table_name__

    def callback(resp: None) -> None:
        pass

    query = prepare_insert_query(table_name=table_name,
                                 columns_names=instance.keys())
    prepared_statement = session.prepare(query)
    future = session.execute_async(prepared_statement,
                                   parameters=instance.values())
    future.add_callback(callback)
    for num, instance in enumerate(instances, start=2):
        future = session.execute_async(prepared_statement,
                                       parameters=instance.values())
        future.add_callback(callback)
    msg = (f'Successfully finished processing '
           f'"{table_name}" table\'s records, '
           f'number of instances processed: {num}.')
    logger.info(msg)


def prepare_insert_query(table_name: str,
                         columns_names: List[str]) -> str:
    binds_aliases = [PREPARED_BIND_ALIAS] * len(columns_names)
    columns_names_str = join_str(columns_names)
    binds_aliases_str = join_str(binds_aliases)
    return (f'INSERT INTO {table_name} ({columns_names_str}) '
            f'VALUES ({binds_aliases_str})')

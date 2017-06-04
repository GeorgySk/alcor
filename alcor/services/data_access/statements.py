from typing import List

from cassandra.cqlengine.models import Model
from cassandra.cqlengine.query import AbstractQuerySet
from cassandra.cqlengine.statements import SelectStatement

from alcor.utils import join_str

ALL_COLUMNS_ALIAS = '*'

BIND_ALIAS = '?'


def select_statement(*,
                     table_name: str,
                     columns_names: List[str] = ALL_COLUMNS_ALIAS,
                     limit: int = None):
    columns_names_str = join_str(columns_names)
    statement = (f'SELECT {columns_names_str} '
                 f'FROM {table_name} ')
    if limit is not None:
        statement += f'LIMIT {limit} '
    return statement


def insert_statement(*,
                     table_name: str,
                     columns_names: List[str]) -> str:
    binds_aliases = [BIND_ALIAS] * len(columns_names)
    columns_names_str = join_str(columns_names)
    binds_aliases_str = join_str(binds_aliases)
    return (f'INSERT INTO {table_name} ({columns_names_str}) '
            f'VALUES ({binds_aliases_str})')


def query_to_select_statement(query: AbstractQuerySet,
                              *,
                              count: bool = False) -> SelectStatement:
    select_statement = query._select_query()
    select_statement.count = count
    return select_statement


def model_insert_statement(model: Model) -> str:
    table_name = model.column_family_name()
    columns_names = model().keys()
    return insert_statement(table_name=table_name,
                            columns_names=columns_names)

import copy
from typing import List, Iterable

import pydevd
from cassandra.cqlengine.functions import QueryValue
from cassandra.cqlengine.models import (Model, ColumnQueryEvaluator)
from cassandra.cqlengine.query import AbstractQuerySet
from cassandra.cqlengine.statements import (SelectStatement, WhereClause)

from alcor.utils import join_str

ASSIGNMENT_CHARACTER = '='

BIND_ALIAS = '?'


def insert_statement(*,
                     table_name: str,
                     columns_names: List[str]) -> str:
    columns_names_str = join_str(columns_names)
    binds_aliases = generate_binds_aliases(columns_names)
    binds_aliases_str = join_str(binds_aliases)
    return (f'INSERT INTO {table_name} ({columns_names_str}) '
            f'VALUES ({binds_aliases_str})')


def update_statement(*,
                     table_name: str,
                     columns_names: List[str],
                     where_clauses: List[WhereClause]) -> str:
    binds_aliases = generate_binds_aliases(columns_names)
    assignments = map(ASSIGNMENT_CHARACTER.join, zip(columns_names,
                                                     binds_aliases))
    assignments_str = join_str(assignments)
    modified_where_clauses = [copy.deepcopy(where_clause)
                              for where_clause in where_clauses]
    for where_clause in modified_where_clauses:
        value = where_clause.value
        query_value = QueryValue(value)
        if isinstance(value, str):
            format_string = '\'{0}\''
        else:
            format_string = '{0}'
        query_value.format_string = format_string
        query_value.context_id = value
        where_clause.query_value = query_value
    where_str = join_str(modified_where_clauses, sep=' AND ')
    return (f'UPDATE {table_name} '
            f'SET {assignments_str} '
            f'WHERE {where_str}')


def generate_binds_aliases(columns_names: List[str]) -> List[str]:
    return [BIND_ALIAS] * len(columns_names)


def query_to_select_statement(query: AbstractQuerySet,
                              *,
                              count: bool = False) -> SelectStatement:
    select_statement = query._select_query()
    select_statement.count = count
    return select_statement


def model_insert_statement(model: Model,
                           *,
                           include_keyspace: bool = True) -> str:
    table_name = model.column_family_name(include_keyspace)
    columns_names = model().keys()
    return insert_statement(table_name=table_name,
                            columns_names=columns_names)


def model_update_statement(model: Model,
                           *,
                           columns: Iterable[ColumnQueryEvaluator] = None,
                           where_clauses: List[WhereClause],
                           include_keyspace: bool = True) -> str:
    table_name = model.column_family_name(include_keyspace)
    if columns:
        columns_names = [column.column.db_field_name
                         for column in columns]
    else:
        columns_names = model().keys()
    return update_statement(table_name=table_name,
                            columns_names=columns_names,
                            where_clauses=where_clauses)

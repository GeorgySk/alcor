from typing import (Iterable,
                    List)

from cassandra.cqlengine.models import (Model,
                                        ColumnQueryEvaluator)
from cassandra.cqlengine.operators import InOperator
from cassandra.cqlengine.query import AbstractQuerySet
from cassandra.cqlengine.statements import (SelectStatement,
                                            WhereClause,
                                            InQuoter)
from cassandra.encoder import cql_quote

from alcor.utils import join_str

ASSIGNMENT_CHARACTER = '='

BIND_ALIAS = '?'


def query_to_select_statement(query: AbstractQuerySet,
                              *,
                              count: bool = False
                              ) -> SelectStatement:
    select_statement = query._select_query()
    select_statement.count = count
    return select_statement


def model_insert_statement(model: Model,
                           *,
                           columns: Iterable[ColumnQueryEvaluator] = None,
                           include_keyspace: bool = True) -> str:
    table_name = model.column_family_name(include_keyspace)
    columns_names = model_columns_names(model,
                                        columns=columns)
    return insert_statement(table_name=table_name,
                            columns_names=columns_names)


def model_update_statement(model: Model,
                           *,
                           columns: Iterable[ColumnQueryEvaluator] = None,
                           where_clauses: List[WhereClause],
                           include_keyspace: bool = True) -> str:
    table_name = model.column_family_name(include_keyspace)
    columns_names = model_columns_names(model,
                                        columns=columns)
    return update_statement(table_name=table_name,
                            columns_names=columns_names,
                            where_clauses=where_clauses)


def model_columns_names(model,
                        *,
                        columns: Iterable[ColumnQueryEvaluator] = None
                        ) -> List[str]:
    # FIXME: find a better way of getting model's columns names
    if columns:
        return [column.column.db_field_name
                for column in columns]
    return list(model._columns.keys())


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
    where_clauses_str = where_clauses_to_str(where_clauses)
    return (f'UPDATE {table_name} '
            f'SET {assignments_str} '
            f'WHERE {where_clauses_str}')


def where_clauses_to_str(clauses: List[WhereClause]
                         ) -> str:
    quoted_values = (
        InQuoter(clause.value)
        if isinstance(clause.operator,
                      InOperator)
        else cql_quote(clause.value)
        for clause in clauses)
    clauses_strs = (
        '{field} {operator} {value}'.format(
            field=clause.field,
            operator=clause.operator,
            value=value)
        for clause, value in zip(clauses,
                                 quoted_values))
    return join_str(clauses_strs,
                    sep=' AND ')


def generate_binds_aliases(columns_names: List[str]) -> List[str]:
    return [BIND_ALIAS] * len(columns_names)

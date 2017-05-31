from typing import List

from alcor.utils import join_str

ALL_COLUMNS_ALIAS = '*'

BIND_ALIAS = '?'


def select_statement(*,
                     table_name: str,
                     columns_names: str = ALL_COLUMNS_ALIAS,
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

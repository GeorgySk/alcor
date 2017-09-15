import logging
import os
from typing import (Union,
                    Callable,
                    Iterator,
                    Dict,
                    Tuple,
                    List)

import numpy as np

from alcor.services.simulations import (da_color,
                                        db_color,
                                        one_table,
                                        da_cooling,
                                        db_cooling)
from alcor.types import CoolingSequencesType

logger = logging.getLogger(__name__)


def read(table_name: str) -> Union[CoolingSequencesType,
                                   Dict[str, np.ndarray]]:
    logging.basicConfig(format='%(filename)s %(funcName)s '
                               '%(levelname)s: %(message)s',
                        level=logging.DEBUG)

    tables = {'da_color': da_color,
              'db_color': db_color,
              'one_table': one_table,
              'da_cooling': da_cooling,
              'db_cooling': db_cooling}

    try:
        table = tables[table_name]
    except KeyError as err:
        err_msg = f'Invalid table name: "{table_name}", not found'
        raise ValueError(err_msg) from err
    split_table_by_metallicities = table.split_by_metallicities

    if split_table_by_metallicities:
        return filled_table_split_by_metallicity(
            folders=table.FILES_FOLDER,
            files_paths_by_metallicities=table.FILES_PATHS,
            metallicities_per_thousand=table.METALLICITIES_PER_THOUSAND,
            masses=table.MASSES,
            pre_wd_lifetimes=table.PRE_WD_LIFETIMES,
            rows_count=table.ROWS_COUNT,
            fill_types_by_metallicities=table.FILL_TYPES_BY_METALLICITIES,
            fill_rules=table.SEQUENCE_FILL_RULES)
    else:
        return filled_table(folder=table.FILES_FOLDER,
                            files_paths=table.FILES_PATHS,
                            rows_count=table.ROWS_COUNT,
                            masses=table.MASSES,
                            fill_rules=table.SEQUENCE_FILL_RULES)


def filled_table(folder: str,
                 files_paths: List[str],
                 rows_count: int,
                 masses: np.ndarray,
                 fill_rules: Dict[str, Dict[str, int]]
                 ) -> Dict[str, np.ndarray]:
    files_paths = [
        os.path.join(folder, file_path)
        for file_path in files_paths]

    files_count = len(files_paths)
    shape = (files_count, rows_count)
    table = dict(mass=masses,
                 rows_counts=np.empty(files_count, dtype='i'))

    for sequence in fill_rules.keys():
        table[sequence] = nan_matrix(shape)

    read_files(files_paths=files_paths,
               table=table,
               max_rows=rows_count,
               fill_rules=fill_rules)

    return table


def filled_table_split_by_metallicity(
        folders: Dict[int, str],
        files_paths_by_metallicities: Dict[int, List[str]],
        metallicities_per_thousand: List[int],
        masses: List[np.ndarray],
        pre_wd_lifetimes: List[np.ndarray],
        rows_count: int,
        fill_types_by_metallicities: Dict[int, int],
        fill_rules: Dict[str, Dict[str, Union[int, List[int], Callable]]]
        ) -> CoolingSequencesType:
    for metallicity, folder_path in folders.items():
        files_paths = files_paths_by_metallicities[metallicity]
        files_paths_by_metallicities[metallicity] = [
            os.path.join(folder_path, file_path)
            for file_path in files_paths]

    table = dict(metallicities_cooling_sequences(
        metallicities_per_thousand,
        files_paths_by_metallicities,
        masses,
        pre_wd_lifetimes,
        rows_count))

    for metallicity, fill_type in (fill_types_by_metallicities.items()):
        read_files(
            files_paths=files_paths_by_metallicities[metallicity],
            table=table[metallicity],
            fill_type=fill_type,
            max_rows=rows_count,
            fill_rules=fill_rules)

    return table


def metallicities_cooling_sequences(
        metallicities: List[int],
        files_paths_by_metallicities: Dict[int, List[str]],
        masses: List[np.ndarray],
        pre_wd_lifetimes: List[np.ndarray],
        rows_count: int
        ) -> Iterator[Tuple[int, Dict[str, np.ndarray]]]:
    for metallicity, mass, pre_wd_lifetime in zip(metallicities,
                                                  masses,
                                                  pre_wd_lifetimes):
        # TODO: I think I don't need files paths. length of masses is enough
        files_count = len(files_paths_by_metallicities[metallicity])
        shape = (files_count, rows_count)
        cooling_sequence = dict(mass=mass,
                                pre_wd_lifetime=pre_wd_lifetime,
                                cooling_time=nan_matrix(shape),
                                effective_temperature=nan_matrix(shape),
                                surface_gravity=nan_matrix(shape),
                                luminosity=nan_matrix(shape),
                                rows_counts=np.empty(files_count, dtype='i'))
        yield metallicity, cooling_sequence


def nan_matrix(shape: Tuple[int, ...]) -> np.ndarray:
    return np.full(shape, np.nan)


def read_files(files_paths: List[str],
               table: Dict[str, np.ndarray],
               max_rows: int,
               fill_rules:
                   Dict[str, Dict[str, Union[int, List[int], Callable]]],
               fill_type: int = None) -> None:
    rows_counts = table['rows_counts']

    for file_path_index, file_path in enumerate(files_paths):
        with open(file_path, 'r') as file:
            lines = [line.split() for line in file]
            for row_index, row in enumerate(lines):
                if row_index == max_rows:
                    break
                # In Fortran indexation starts from 1
                rows_counts[file_path_index] = row_index + 1
                for sequence, fill_rule in fill_rules.items():
                    if fill_type is None:
                        if 'converting_method' in fill_rule:
                            table[sequence][file_path_index, row_index] = (
                                fill_rule['converting_method'](
                                    float(row[fill_rule['column']])))
                        else:
                            table[sequence][file_path_index, row_index] = (
                                float(row[fill_rule['column']]))
                    else:
                        if 'converting_method' in fill_rule:
                            table[sequence][file_path_index, row_index] = (
                                fill_rule['converting_method'](
                                    float(row[fill_rule['column'][fill_type]]),
                                    fill_type=fill_type,
                                    pre_wd_lifetime=table['pre_wd_lifetime'][
                                        file_path_index]))
                        else:
                            table[sequence][file_path_index, row_index] = (
                                float(row[fill_rule['column'][fill_type]]))

if __name__ == '__main__':
    read(table_name='da_cooling')

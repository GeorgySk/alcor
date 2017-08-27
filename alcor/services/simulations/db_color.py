import csv
import os
from typing import (Dict,
                    Tuple,
                    List)

import numpy as np


def initialize_sequences(rows_count: int = 60) -> None:
    files_paths = ['Table_Mass_05_DB_sort.trk',
                   'Table_Mass_06_DB_sort.trk',
                   'Table_Mass_07_DB_sort.trk',
                   'Table_Mass_08_DB_sort.trk',
                   'Table_Mass_09_DB_sort.trk',
                   'Table_Mass_10_DB_sort.trk',
                   'Table_Mass_12_DB_sort.trk']
    files_folder = './db_color_tables'

    base_dir = os.path.dirname(__file__)

    files_paths = [
        os.path.join(base_dir, files_folder, file_path)
        for file_path in files_paths]

    masses = [np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2])]

    files_count = len(files_paths)
    shape = (files_count, rows_count)
    color_table = dict(mass=masses,
                       luminosity=nan_matrix(shape),
                       u_ubvri_absolute=nan_matrix(shape),
                       b_ubvri_absolute=nan_matrix(shape),
                       v_ubvri_absolute=nan_matrix(shape),
                       r_ubvri_absolute=nan_matrix(shape),
                       i_ubvri_absolute=nan_matrix(shape),
                       rows_counts=np.empty(files_count, dtype='i'))
    read_files(files_paths=files_paths,
               color_table=color_table,
               max_rows=rows_count)


def nan_matrix(shape: Tuple[int, ...]) -> np.ndarray:
    return np.full(shape, np.nan)


def read_files(files_paths: List[str],
               color_table: Dict[str, np.ndarray],
               max_rows: int) -> None:
    luminosity = color_table['luminosity']
    u_ubvri_absolute = color_table['u_ubvri_absolute']
    b_ubvri_absolute = color_table['b_ubvri_absolute']
    v_ubvri_absolute = color_table['v_ubvri_absolute']
    r_ubvri_absolute = color_table['r_ubvri_absolute']
    i_ubvri_absolute = color_table['i_ubvri_absolute']
    rows_counts = color_table['rows_counts']

    for file_path_index, file_path in enumerate(files_paths):
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file,
                                    delimiter=' ',
                                    skipinitialspace=True)
            for row_index, row in enumerate(csv_reader):
                if row_index == max_rows:
                    break
                # In Fortran indexation starts from 1
                rows_counts[file_path_index] = row_index + 1
                luminosity[file_path_index, row_index] = -(float(row[2])
                                                           - 4.72) / 2.5
                u_ubvri_absolute[file_path_index, row_index] = float(row[4])
                b_ubvri_absolute[file_path_index, row_index] = float(row[5])
                v_ubvri_absolute[file_path_index, row_index] = float(row[6])
                r_ubvri_absolute[file_path_index, row_index] = float(row[7])
                i_ubvri_absolute[file_path_index, row_index] = float(row[8])

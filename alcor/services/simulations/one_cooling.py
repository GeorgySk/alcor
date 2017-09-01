import csv
import os
from math import log10
from typing import (Dict,
                    Tuple,
                    List)

import numpy as np


# TODO: give better names to everything
def initialize_sequences(rows_count: int = 300
                         ) -> None:
    color_tables_files_paths = ['color_106.out',
                                'color_110.out',
                                'color_116.out',
                                'color_120.out',
                                'color_124.out',
                                'color_128.out']
    cooling_tables_files_paths = ['t106_he.trk',
                                  't110_he.trk',
                                  't116_he.trk',
                                  't120_he.trk',
                                  't128_he.trk']
    # TODO: maybe I should place everything in one folder
    cooling_tables_folder = './one_color_tables'
    color_tables_folder = './one_cooling_tables'

    base_dir = os.path.dirname(__file__)

    color_tables_files_paths = [
        os.path.join(base_dir, color_tables_folder, file_path)
        for file_path in color_tables_files_paths]

    cooling_tables_files_paths = [
        os.path.join(base_dir, cooling_tables_folder, file_path)
        for file_path in cooling_tables_files_paths]

    color_masses = [np.array([1.06, 1.10, 1.16, 1.20, 1.24, 1.28])]
    color_pre_wd_lifetimes = [np.zeros(len(color_tables_files_paths))]
    cooling_masses = [np.array([1.06, 1.10, 1.16, 1.20, 1.28])]
    cooling_pre_wd_lifetimes = [np.zeros(len(cooling_tables_files_paths))]

    files_count = len(color_tables_files_paths)
    shape = (files_count, rows_count)
    color_table = dict(mass=color_masses,
                       pre_wd_lifetimes=color_pre_wd_lifetimes,
                       luminosity=nan_matrix(shape),
                       v_ubvri_absolute=nan_matrix(shape),
                       log_cooling_time=nan_matrix(shape),
                       log_effective_temperature=nan_matrix(shape),
                       bv_ubvri=nan_matrix(shape),
                       vi_ubvri=nan_matrix(shape),
                       vr_ubvri=nan_matrix(shape),
                       uv_ubvri=nan_matrix(shape),
                       rows_counts=np.empty(files_count, dtype='i'))
    read_color_files(
        files_paths=color_tables_files_paths,
        cooling_sequence=color_table,
        max_rows=rows_count)

    files_count = len(cooling_tables_files_paths)
    shape = (files_count, rows_count)
    cooling_table = dict(mass=cooling_masses,
                         pre_wd_lifetimes=cooling_pre_wd_lifetimes,
                         log_surface_gravity=nan_matrix(shape),
                         log_cooling_time=nan_matrix(shape),
                         rows_counts=np.empty(files_count, dtype='i'))
    read_cooling_files(
        files_paths=cooling_tables_files_paths,
        cooling_sequence=cooling_table,
        max_rows=rows_count)


def nan_matrix(shape: Tuple[int, ...]) -> np.ndarray:
    return np.full(shape, np.nan)


def read_color_files(files_paths: List[str],
                     cooling_sequence: Dict[str, np.ndarray],
                     max_rows: int) -> None:
    luminosity = cooling_sequence['luminosity']
    v_ubvri_absolute = cooling_sequence['v_ubvri_absolute']
    log_cooling_time = cooling_sequence['log_cooling_time']
    log_effective_temperature = cooling_sequence['log_effective_temperature']
    bv_ubvri = cooling_sequence['bv_ubvri']
    vi_ubvri = cooling_sequence['vi_ubvri']
    vr_ubvri = cooling_sequence['vr_ubvri']
    uv_ubvri = cooling_sequence['uv_ubvri']
    rows_counts = cooling_sequence['rows_counts']

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
                luminosity[file_path_index, row_index] = float(row[0])
                v_ubvri_absolute[file_path_index, row_index] = float(row[11])
                log_cooling_time[file_path_index, row_index] = float(row[12])
                log_effective_temperature[file_path_index, row_index] = float(
                    row[1])
                bv_ubvri[file_path_index, row_index] = float(row[2])
                vi_ubvri[file_path_index, row_index] = float(row[8])
                vr_ubvri[file_path_index, row_index] = float(row[3])
                uv_ubvri[file_path_index, row_index] = float(row[9])


def read_cooling_files(files_paths: List[str],
                       cooling_sequence: Dict[str, np.ndarray],
                       max_rows: int) -> None:
    log_surface_gravity = cooling_sequence['log_surface_gravity']
    log_cooling_time = cooling_sequence['log_cooling_time ']
    rows_counts = cooling_sequence['rows_counts']

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
                log_surface_gravity[file_path_index, row_index] = log10(
                    10. ** float(row[2]) / 6.96E10)
                log_cooling_time[file_path_index, row_index] = float(row[4])

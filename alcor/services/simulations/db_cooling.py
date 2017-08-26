import csv
import os
from typing import (Iterator,
                    Dict,
                    Tuple,
                    List)

import numpy as np

from alcor.types import CoolingSequencesType


def initialize_sequences(rows_count: int = 400
                         ) -> CoolingSequencesType:
    # Metallicities were multiplied by 1000 in order to keep dict.keys as ints
    files_paths_by_metallicities = {1: ['05047_db_Z=0.001',
                                        '05527_db_Z=0.001',
                                        '059328_db_Z=0.001',
                                        '062738_db_Z=0.001',
                                        '06602_db_Z=0.001',
                                        '069289_db_Z=0.001',
                                        '08637_db_Z=0.001'],
                                    10: ['db_cool_0514.seq',
                                         'db_cool_0530.seq',
                                         'db_cool_0542.seq',
                                         'db_cool_0565.seq',
                                         'db_cool_0584.seq',
                                         'db_cool_0609.seq',
                                         'db_cool_0664.seq',
                                         'db_cool_0741.seq',
                                         'db_cool_0869.seq'],
                                    60: ['0524db_006_sflhdiff.trk',
                                         '0570db_006_sflhdiff.trk',
                                         '0593db_006_sflhdiff.trk',
                                         '061db_006_sflhdiff.trk',
                                         '0632db_006_sflhdiff.trk',
                                         '0659db_006_sflhdiff.trk',
                                         '070db_006_sflhdiff.trk',
                                         '076db_006_sflhdiff.trk',
                                         '087db_006_sflhdiff.trk']}
    folders_paths_by_metallicities = {
        1:  './db_cooling_tables/Z0001',
        10: './db_cooling_tables/Z001',
        60: './db_cooling_tables/Z006'}

    base_dir = os.path.dirname(__file__)

    for metallicity, folder_path in folders_paths_by_metallicities.items():
        files_paths = files_paths_by_metallicities[metallicity]
        files_paths_by_metallicities[metallicity] = [
            os.path.join(base_dir, folder_path, file_path)
            for file_path in files_paths]

    metallicities_per_thousand = [1, 10, 60]
    masses = [np.array([0.5047, 0.5527, 0.59328, 0.62738,
                        0.6602, 0.69289, 0.8637]),
              np.array([0.514, 0.53, 0.542, 0.565, 0.584,
                        0.61, 0.664, 0.741, 0.869]),
              np.array([0.524, 0.570, 0.593, 0.61, 0.632,
                        0.659, 0.70, 0.76, 0.87])]
    pre_wd_lifetimes = [
        np.zeros(len(files_paths_by_metallicities[1])),
        np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892,
                  0.7422, 0.4431, 0.287, 0.114]),
        np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892,
                  0.7422, 0.4431, 0.287, 0.114])]

    cooling_sequences_by_metallicities = dict(metallicities_cooling_sequences(
        metallicities_per_thousand,
        files_paths_by_metallicities,
        masses,
        pre_wd_lifetimes,
        rows_count))

    fill_types_by_metallicities = {1: 1,
                                   10: 2,
                                   60: 1}

    for metallicity, fill_type in fill_types_by_metallicities.items():
        read_files(
            files_paths=files_paths_by_metallicities[metallicity],
            cooling_sequence=cooling_sequences_by_metallicities[metallicity],
            fill_type=fill_type,
            max_rows=rows_count)

    return cooling_sequences_by_metallicities


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
               cooling_sequence: Dict[str, np.ndarray],
               fill_type: int,
               max_rows: int) -> None:
    cooling_time = cooling_sequence['cooling_time']
    effective_temperature = cooling_sequence['effective_temperature']
    surface_gravity = cooling_sequence['surface_gravity']
    luminosity = cooling_sequence['luminosity']
    pre_wd_lifetime = cooling_sequence['pre_wd_lifetime']
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
                effective_temperature[file_path_index, row_index] = (
                    10. ** float(row[1]))
                if fill_type == 1:
                    cooling_time[file_path_index, row_index] = (
                        10. ** float(row[8]) / 1000.
                        - pre_wd_lifetime[file_path_index])
                    surface_gravity[file_path_index, row_index] = float(
                        row[22])
                if fill_type == 2:
                    cooling_time[file_path_index, row_index] = float(row[2])
                    surface_gravity[file_path_index, row_index] = float(row[3])

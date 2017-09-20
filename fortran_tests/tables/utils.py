from math import isclose
from typing import (Union,
                    Callable,
                    Tuple,
                    List)

import numpy as np

from alcor.types import CoolingSequencesType


# TODO: rename rows_count to max_rows_count
def read_sequences_by_metallicity_from_fortran(
        rows_count: int,
        files_counts_per_metallicity: List[int],
        fill_types: List[int],
        fort_files_initial_units: List[int],
        metallicities_by_thousand: List[int],
        get_from_fortran: Callable,
        ) -> CoolingSequencesType:
    metallicities_count = len(files_counts_per_metallicity)
    rows_counts = [nan_matrix(shape=files_count, dtype='i')
                   for files_count in files_counts_per_metallicity]
    masses = [nan_matrix(shape=files_count)
              for files_count in files_counts_per_metallicity]
    cooling_times = [nan_matrix(shape=(files_count, rows_count), order='F')
                     for files_count in files_counts_per_metallicity]
    pre_wd_lifetimes = [nan_matrix(shape=files_count)
                        for files_count in files_counts_per_metallicity]
    luminosities = [nan_matrix(shape=(files_count, rows_count), order='F')
                    for files_count in files_counts_per_metallicity]
    effective_temperatures = [
        nan_matrix(shape=(files_count, rows_count), order='F')
        for files_count in files_counts_per_metallicity]
    surface_gravities = [
        nan_matrix(shape=(files_count, rows_count), order='F')
        for files_count in files_counts_per_metallicity]

    for index in range(metallicities_count):
        get_from_fortran(flag=fill_types[index],
                         initlink=fort_files_initial_units[index],
                         ntrk=rows_counts[index],
                         mass=masses[index],
                         coolingtime=cooling_times[index],
                         prevtime=pre_wd_lifetimes[index],
                         luminosity=luminosities[index],
                         efftemp=effective_temperatures[index],
                         gravacc=surface_gravities[index])

    cooling_sequences_by_metallicities = {
        metallicities_by_thousand[index]: dict(
            mass=masses[index],
            pre_wd_lifetime=pre_wd_lifetimes[index],
            cooling_time=cooling_times[index],
            effective_temperature=effective_temperatures[index],
            surface_gravity=surface_gravities[index],
            luminosity=luminosities[index],
            rows_counts=rows_counts[index])
        for index in range(metallicities_count)}

    return cooling_sequences_by_metallicities


def nan_matrix(shape: Union[int, Tuple[int, ...]],
               dtype: str = 'f',
               order: str = 'C'
               ) -> np.ndarray:
    return np.full(shape=shape,
                   fill_value=np.nan,
                   dtype=dtype,
                   order=order)


def values_by_metallicity_are_close(table: CoolingSequencesType,
                                    other_table: CoolingSequencesType,
                                    *,
                                    relative_tolerance: float = 1E-4) -> bool:
    for metallicity, table_sequences in table.items():
        other_table_sequences = other_table[metallicity]
        for sequence_name, sequence_values in table_sequences.items():
            other_sequence_values = other_table_sequences[sequence_name]

            if sequence_values.ndim == 1:
                for x_value, y_value in zip(sequence_values,
                                            other_sequence_values):
                    # NaN when compared with other NaN always gives False
                    if np.isnan(x_value) and np.isnan(y_value):
                        continue
                    if not isclose(x_value,
                                   y_value,
                                   rel_tol=relative_tolerance):
                        return False
            else:
                for x_row, y_row in zip(sequence_values,
                                        other_sequence_values):
                    for x_value, y_value in zip(x_row, y_row):
                        if np.isnan(x_value) and np.isnan(y_value):
                            continue
                        if not isclose(x_value,
                                       y_value,
                                       rel_tol=relative_tolerance):
                            return False
    return True

import logging
import operator
from functools import reduce
from math import isclose
from typing import (Any,
                    Union,
                    Callable,
                    Mapping,
                    Hashable,
                    Tuple,
                    List)

import numpy as np

from alcor.types import CoolingSequencesType

logger = logging.getLogger(__name__)


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
    rows_counts = [nan_matrix(shape=files_count,
                              dtype='i')
                   for files_count in files_counts_per_metallicity]
    masses = [nan_matrix(shape=files_count)
              for files_count in files_counts_per_metallicity]
    cooling_times = [nan_matrix(shape=(files_count, rows_count),
                                order='F')
                     for files_count in files_counts_per_metallicity]
    pre_wd_lifetimes = [nan_matrix(shape=files_count)
                        for files_count in files_counts_per_metallicity]
    luminosities = [nan_matrix(shape=(files_count, rows_count),
                               order='F')
                    for files_count in files_counts_per_metallicity]
    effective_temperatures = [
        nan_matrix(shape=(files_count, rows_count),
                   order='F')
        for files_count in files_counts_per_metallicity]
    surface_gravities = [
        nan_matrix(shape=(files_count, rows_count),
                   order='F')
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
    for _, (table_sequences,
            other_table_sequences) in zip_mappings(table,
                                                   other_table):
        for sequence_name, (sequence_values,
                            other_sequence_values) in zip_mappings(
                table_sequences,
                other_table_sequences):
            if sequence_values.ndim == 1 and other_sequence_values.ndim == 1:
                if not one_dim_arrays_are_close(
                        sequence_values,
                        other_sequence_values,
                        relative_tolerance=relative_tolerance):
                    return False
            elif sequence_values.ndim == 2 and other_sequence_values.ndim == 2:
                if not two_dim_arrays_are_close(
                        sequence_values,
                        other_sequence_values,
                        relative_tolerance=relative_tolerance):
                    return False
            else:
                logger.error(f'Dimensions mismatch for {sequence_name}')
                return False
    return True


def zip_mappings(*mappings: Mapping[Hashable, Any]):
    keys_sets = map(set, mappings)
    common_keys = reduce(set.intersection, keys_sets)
    for key in common_keys:
        yield key, tuple(map(operator.itemgetter(key), mappings))


# TODO: check if it's achievable by numpy functions
def one_dim_arrays_are_close(x_values: np.ndarray,
                             y_values: np.ndarray,
                             *,
                             relative_tolerance: float) -> bool:
    # TODO: check for zip alternatives
    for x_value, y_value in zip(x_values,
                                y_values):
        # NaN when compared with other NaN always gives False
        if np.isnan(x_value) and np.isnan(y_value):
            continue
        if not isclose(x_value,
                       y_value,
                       rel_tol=relative_tolerance):
            return False
    return True


# TODO: is it possible to use np.reduce here?
def two_dim_arrays_are_close(x_matrix: np.ndarray,
                             y_matrix: np.ndarray,
                             *,
                             relative_tolerance: float) -> bool:
    for x_row, y_row in zip(x_matrix,
                            y_matrix):
        if not one_dim_arrays_are_close(x_row,
                                        y_row,
                                        relative_tolerance=relative_tolerance):
            return False
    return True

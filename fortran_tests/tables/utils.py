import logging
from math import isclose
from typing import (Union,
                    Callable,
                    Dict,
                    Tuple,
                    List)

import numpy as np

from alcor.types import CoolingSequencesType
from alcor.utils import zip_mappings

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


# TODO: rename rows_count to max_rows_count
def read_colors_from_fortran(rows_count: int,
                             files_count: int,
                             fort_files_initial_unit: int,
                             get_from_fortran: Callable
                             ) -> Dict[str, np.ndarray]:
    rows_counts = nan_matrix(shape=files_count, dtype='i')
    masses = nan_matrix(shape=files_count)
    luminosities = nan_matrix(shape=(files_count, rows_count),
                              order='F')
    u_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')
    b_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')
    v_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')
    r_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')
    i_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')
    j_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')

    get_from_fortran(initlink=fort_files_initial_unit,
                     ntrk=rows_counts,
                     mass=masses,
                     luminosity=luminosities,
                     color_u=u_ubvri,
                     color_b=b_ubvri,
                     color_v=v_ubvri,
                     color_r=r_ubvri,
                     color_i=i_ubvri,
                     color_j=j_ubvri)

    return dict(mass=masses,
                luminosity=luminosities,
                u_ubvri_absolute=u_ubvri,
                b_ubvri_absolute=b_ubvri,
                v_ubvri_absolute=v_ubvri,
                r_ubvri_absolute=r_ubvri,
                i_ubvri_absolute=i_ubvri,
                j_ubvri_absolute=j_ubvri,
                rows_counts=rows_counts)


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


def values_are_close(x: Dict[str, np.ndarray],
                     y: Dict[str, np.ndarray],
                     rel_tol: float = 1E-4) -> bool:
    for key, x_values in x.items():
        y_values = y[key]
        for x_value, y_value in zip(x_values, y_values):
            if x_value.shape is ():
                if not np.isnan(x_value) and not np.isnan(y_value):
                    if not isclose(x_value,
                                   y_value,
                                   rel_tol=rel_tol):
                        return False
            else:
                for (x_sub_value, y_sub_value) in zip(x_value, y_value):
                    if (not np.isnan(x_sub_value)
                            and not np.isnan(y_sub_value)):
                        if not isclose(x_sub_value,
                                       y_sub_value,
                                       rel_tol=rel_tol):
                            return False
    return True


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

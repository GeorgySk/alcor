from math import isclose
from typing import (Union,
                    Callable,
                    Dict,
                    Tuple,
                    List)

from alcor.services.simulations import table
from alcor.types import CoolingSequencesType

import numpy as np

import da_cooling
import db_cooling
import da_color


def test_table() -> None:
    da_cooling_table_by_python = table.read(table_name='da_cooling')
    da_cooling_table_by_fortran = read_sequences_by_metallicity_from_fortran(
        rows_count=650,
        files_counts_per_metallicity=[7, 10, 8, 8],
        fill_types=[1, 2, 3, 3],
        fort_files_initial_units=[10, 20, 30, 40],
        metallicities_by_thousand=[1, 10, 30, 60],
        get_from_fortran=da_cooling.incoolda)

    assert values_by_metallicity_are_close(da_cooling_table_by_python,
                                           da_cooling_table_by_fortran)

    db_cooling_table_by_python = table.read(table_name='db_cooling')
    db_cooling_table_by_fortran = read_sequences_by_metallicity_from_fortran(
        rows_count=400,
        files_counts_per_metallicity=[7, 9, 9],
        fill_types=[1, 2, 3],
        fort_files_initial_units=[90, 100, 110],
        metallicities_by_thousand=[1, 10, 60],
        get_from_fortran=db_cooling.incooldb)

    assert values_by_metallicity_are_close(db_cooling_table_by_python,
                                           db_cooling_table_by_fortran)

    da_color_table_by_python = table.read(table_name='da_color')
    da_color_table_by_fortran = read_colors_from_fortran(
        rows_count=650,
        files_count=10,
        fort_files_initial_unit=60,
        get_from_fortran=da_color.color)

    assert values_are_close(da_color_table_by_python,
                            da_color_table_by_fortran)


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


# TODO: rename rows_count to max_rows_count
def read_colors_from_fortran(rows_count: int,
                             files_count: int,
                             fort_files_initial_unit: int,
                             get_from_fortran: Callable
                             ) -> Dict[str, np.ndarray]:
    rows_counts = nan_matrix(shape=files_count, dtype='i')
    masses = nan_matrix(shape=files_count)
    luminosities = nan_matrix(shape=(files_count, rows_count), order='F')
    u_ubvri = nan_matrix(shape=(files_count, rows_count), order='F')
    b_ubvri = nan_matrix(shape=(files_count, rows_count), order='F')
    v_ubvri = nan_matrix(shape=(files_count, rows_count), order='F')
    r_ubvri = nan_matrix(shape=(files_count, rows_count), order='F')
    i_ubvri = nan_matrix(shape=(files_count, rows_count), order='F')

    get_from_fortran(initlink=fort_files_initial_unit,
                     ntrk=rows_counts,
                     mass=masses,
                     luminosity=luminosities,
                     color_u=u_ubvri,
                     color_b=b_ubvri,
                     color_v=v_ubvri,
                     color_r=r_ubvri,
                     color_i=i_ubvri)

    return dict(mass=masses,
                luminosity=luminosities,
                u_ubvri_absolute=u_ubvri,
                b_ubvri_absolute=b_ubvri,
                v_ubvri_absolute=v_ubvri,
                r_ubvri_absolute=r_ubvri,
                i_ubvri_absolute=i_ubvri,
                rows_counts=rows_counts)


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


def values_by_metallicity_are_close(x: CoolingSequencesType,
                                    y: CoolingSequencesType,
                                    rel_tol: float = 1E-4) -> bool:
    for metallicity, x_sequences in x.items():
        y_sequences = y[metallicity]
        for key, x_values in x_sequences.items():
            y_values = y_sequences[key]
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

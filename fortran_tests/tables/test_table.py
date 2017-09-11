from math import isclose
from typing import (Union,
                    Tuple)

from alcor.services.simulations import table
from alcor.types import CoolingSequencesType

import numpy as np

import da_cooling
import db_cooling


def test_table() -> None:
    da_cooling_table_by_python = table.read(table_name='da_cooling')
    da_cooling_table_by_fortran = read_sequences_from_fortran()

    assert values_are_close(da_cooling_table_by_python,
                            da_cooling_table_by_fortran)

    db_cooling_table_by_python = table.read(table_name='db_cooling')
    db_cooling_table_by_fortran = read_db_sequences_from_fortran()

    assert values_are_close(db_cooling_table_by_python,
                            db_cooling_table_by_fortran)


def read_db_sequences_from_fortran(rows_count: int = 400
                                   ) -> CoolingSequencesType:
    files_counts_per_metallicity = [7, 9, 9]
    fill_types = [1, 2, 3]
    fort_files_initial_units = [90, 100, 110]
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

    for index in range(3):
        db_cooling.incooldb(flag=fill_types[index],
                            initlink=fort_files_initial_units[index],
                            ntrk=rows_counts[index],
                            mass=masses[index],
                            coolingtime=cooling_times[index],
                            prevtime=pre_wd_lifetimes[index],
                            luminosity=luminosities[index],
                            efftemp=effective_temperatures[index],
                            gravacc=surface_gravities[index])

    metallicities_by_thousand = [1, 10, 60]

    cooling_sequences_by_metallicities = {
        metallicities_by_thousand[index]: dict(
            mass=masses[index],
            pre_wd_lifetime=pre_wd_lifetimes[index],
            cooling_time=cooling_times[index],
            effective_temperature=effective_temperatures[index],
            surface_gravity=surface_gravities[index],
            luminosity=luminosities[index],
            rows_counts=rows_counts[index])
        for index in range(3)}

    return cooling_sequences_by_metallicities


def read_sequences_from_fortran(rows_count: int = 650
                                ) -> CoolingSequencesType:
    files_counts_per_metallicity = [7, 10, 8, 8]
    fill_types = [1, 2, 3, 3]
    fort_files_initial_units = [10, 20, 30, 40]
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

    for index in range(4):
        da_cooling.incoolda(flag=fill_types[index],
                            initlink=fort_files_initial_units[index],
                            ntrk=rows_counts[index],
                            mass=masses[index],
                            coolingtime=cooling_times[index],
                            prevtime=pre_wd_lifetimes[index],
                            luminosity=luminosities[index],
                            efftemp=effective_temperatures[index],
                            gravacc=surface_gravities[index])

    metallicities_by_thousand = [1, 10, 30, 60]

    cooling_sequences_by_metallicities = {
        metallicities_by_thousand[index]: dict(
            mass=masses[index],
            pre_wd_lifetime=pre_wd_lifetimes[index],
            cooling_time=cooling_times[index],
            effective_temperature=effective_temperatures[index],
            surface_gravity=surface_gravities[index],
            luminosity=luminosities[index],
            rows_counts=rows_counts[index])
        for index in range(4)}

    return cooling_sequences_by_metallicities


def nan_matrix(shape: Union[int, Tuple[int, ...]],
               dtype: str = 'f',
               order: str = 'C'
               ) -> np.ndarray:
    return np.full(shape=shape,
                   fill_value=np.nan,
                   dtype=dtype,
                   order=order)


def values_are_close(x: CoolingSequencesType,
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

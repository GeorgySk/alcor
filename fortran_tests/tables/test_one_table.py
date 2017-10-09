from typing import Dict, Callable

import numpy as np

from alcor.services.simulations import table
from alcor.services.simulations.one_table import (ROWS_COUNT,
                                                  FILES_PATHS)
from .utils import (nan_matrix,
                    values_are_close)

import one_tables


def test_one_table() -> None:
    one_table_by_python = table.read(table_name='one_table')
    one_table_by_fortran = read_one_table_from_fortran(
        rows_count=ROWS_COUNT,
        files_count=len(FILES_PATHS),
        get_from_fortran=one_tables.incoolone)

    assert values_are_close(one_table_by_python,
                            one_table_by_fortran)


def read_one_table_from_fortran(rows_count: int,
                                files_count: int,
                                get_from_fortran: Callable
                                ) -> Dict[str, np.ndarray]:
    log_cooling_times = nan_matrix(shape=(files_count, rows_count),
                                   order='F')
    luminosities = nan_matrix(shape=(files_count, rows_count),
                              order='F')
    v_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')
    log_effective_temperatures = nan_matrix(shape=(files_count, rows_count),
                                            order='F')
    masses = nan_matrix(shape=files_count)
    rows_counts = nan_matrix(shape=files_count, dtype='i')
    bv_ubvri = nan_matrix(shape=(files_count, rows_count),
                          order='F')
    vi_ubvri = nan_matrix(shape=(files_count, rows_count),
                          order='F')
    vr_ubvri = nan_matrix(shape=(files_count, rows_count),
                          order='F')
    uv_ubvri = nan_matrix(shape=(files_count, rows_count),
                          order='F')
    j_ubvri = nan_matrix(shape=(files_count, rows_count),
                         order='F')

    get_from_fortran(lgtabone=log_cooling_times,
                     ltabone=luminosities,
                     mvtabone=v_ubvri,
                     lgtetabone=log_effective_temperatures,
                     mtabone=masses,
                     ndatsone=rows_counts,
                     bvtabone=bv_ubvri,
                     vitabone=vi_ubvri,
                     vrtabone=vr_ubvri,
                     uvtabone=uv_ubvri,
                     jone=j_ubvri)

    return dict(mass=masses,
                luminosity=luminosities,
                v_ubvri_absolute=v_ubvri,
                log_cooling_time=log_cooling_times,
                log_effective_temperature=log_effective_temperatures,
                bv_ubvri=bv_ubvri,
                vi_ubvri=vi_ubvri,
                vr_ubvri=vr_ubvri,
                uv_ubvri=uv_ubvri,
                j_ubvri=j_ubvri,
                rows_counts=rows_counts)

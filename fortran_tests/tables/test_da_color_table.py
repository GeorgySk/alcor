from math import isclose
from typing import (Union,
                    Callable,
                    Dict,
                    Tuple)

from alcor.services.simulations import table

import numpy as np

import da_color


def test_da_cooling_table() -> None:
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


def nan_matrix(shape: Union[int, Tuple[int, ...]],
               dtype: str = 'f',
               order: str = 'C'
               ) -> np.ndarray:
    return np.full(shape=shape,
                   fill_value=np.nan,
                   dtype=dtype,
                   order=order)

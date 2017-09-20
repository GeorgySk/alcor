from alcor.services.simulations import table
from .utils import (read_sequences_by_metallicity_from_fortran,
                    values_by_metallicity_are_close)

import da_cooling


def test_da_cooling_table() -> None:
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

from alcor.services.simulations import table
from alcor.services.simulations.da_cooling import (ROWS_COUNT,
                                                   METALLICITIES_PER_THOUSAND,
                                                   FILES_PATHS)
from .utils import (read_sequences_by_metallicity_from_fortran,
                    values_by_metallicity_are_close)

import da_cooling


def test_da_cooling_table() -> None:
    da_cooling_table_by_python = table.read(table_name='da_cooling')

    files_paths = [FILES_PATHS[metallicity]
                   for metallicity in METALLICITIES_PER_THOUSAND]
    files_counts_per_metallicity = list(map(len, files_paths))
    da_cooling_table_by_fortran = read_sequences_by_metallicity_from_fortran(
        rows_count=ROWS_COUNT,
        files_counts_per_metallicity=files_counts_per_metallicity,
        fill_types=[1, 2, 3, 3],
        fort_files_initial_units=[10, 20, 30, 40],
        metallicities_by_thousand=METALLICITIES_PER_THOUSAND,
        get_from_fortran=da_cooling.incoolda)

    assert values_by_metallicity_are_close(da_cooling_table_by_python,
                                           da_cooling_table_by_fortran)

from alcor.services.simulations import table
from alcor.services.simulations.db_cooling import (ROWS_COUNT,
                                                   METALLICITIES_PER_THOUSAND)
from .utils import (read_sequences_by_metallicity_from_fortran,
                    values_by_metallicity_are_close)

import db_cooling


def test_db_cooling_table() -> None:
    db_cooling_table_by_python = table.read(table_name='db_cooling')
    db_cooling_table_by_fortran = read_sequences_by_metallicity_from_fortran(
        rows_count=ROWS_COUNT,
        files_counts_per_metallicity=[7, 9, 9],
        fill_types=[1, 2, 3],
        fort_files_initial_units=[90, 100, 110],
        metallicities_by_thousand=METALLICITIES_PER_THOUSAND,
        get_from_fortran=db_cooling.incooldb)

    assert values_by_metallicity_are_close(db_cooling_table_by_python,
                                           db_cooling_table_by_fortran)

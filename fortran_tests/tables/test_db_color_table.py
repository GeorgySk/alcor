from alcor.services.simulations import table

import db_color

from alcor.services.simulations.db_color import (ROWS_COUNT,
                                                 FILES_PATHS)
from .utils import (read_colors_from_fortran,
                    values_are_close)


def test_db_color_table() -> None:
    db_color_table_by_python = table.read(table_name='db_color')
    db_color_table_by_fortran = read_colors_from_fortran(
        rows_count=ROWS_COUNT,
        files_count=len(FILES_PATHS),
        fort_files_initial_unit=131,
        get_from_fortran=db_color.colordb)

    assert values_are_close(db_color_table_by_python,
                            db_color_table_by_fortran)

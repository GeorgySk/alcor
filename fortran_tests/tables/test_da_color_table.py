from alcor.services.simulations import table

import da_color

from alcor.services.simulations.da_color import (ROWS_COUNT,
                                                 FILES_PATHS)
from .utils import (read_colors_from_fortran,
                    values_are_close)


def test_da_color_table() -> None:
    da_color_table_by_python = table.read(table_name='da_color')
    da_color_table_by_fortran = read_colors_from_fortran(
        rows_count=ROWS_COUNT,
        files_count=len(FILES_PATHS),
        fort_files_initial_unit=60,
        get_from_fortran=da_color.color)

    assert values_are_close(da_color_table_by_python,
                            da_color_table_by_fortran)

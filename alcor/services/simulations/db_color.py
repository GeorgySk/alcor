import numpy as np


def convert_luminosity(luminosity: str) -> float:
    return -(float(luminosity) - 4.72) / 2.5

split_by_metallicities = False

ROWS_COUNT = 60
FILES_PATHS = ['Table_Mass_05_DB_sort.trk',
               'Table_Mass_06_DB_sort.trk',
               'Table_Mass_07_DB_sort.trk',
               'Table_Mass_08_DB_sort.trk',
               'Table_Mass_09_DB_sort.trk',
               'Table_Mass_10_DB_sort.trk',
               'Table_Mass_12_DB_sort.trk']
FILES_DIR = 'db_color_tables'
MASSES = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2])
SEQUENCE_FILL_RULES = {'luminosity': {'column': 2,
                       'converting_method': convert_luminosity},
                       'u_ubvri_absolute': {'column': 4},
                       'b_ubvri_absolute': {'column': 5},
                       'v_ubvri_absolute': {'column': 6},
                       'r_ubvri_absolute': {'column': 7},
                       'i_ubvri_absolute': {'column': 8},
                       'j_ubvri_absolute': {'column': 9}}

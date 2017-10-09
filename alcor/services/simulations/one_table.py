from typing import List

import numpy as np


def convert_j_color(inputs: List[float]) -> float:
    return inputs[0] + inputs[1] - inputs[2] + inputs[3]


split_by_metallicities = False

ROWS_COUNT = 300
FILES_PATHS = ['color_106.out',
               'color_110.out',
               'color_116.out',
               'color_120.out',
               'color_124.out',
               'color_128.out']
FILES_DIR = 'one_tables'
MASSES = np.array([1.06, 1.10, 1.16, 1.20, 1.24, 1.28])
SEQUENCE_FILL_RULES = {'luminosity': {'column': 0},
                       'v_ubvri_absolute': {'column': 11},
                       'log_cooling_time': {'column': 12},
                       'log_effective_temperature': {'column': 1},
                       'bv_ubvri': {'column': 2},
                       'vi_ubvri': {'column': 8},
                       'vr_ubvri': {'column': 3},
                       'uv_ubvri': {'column': 9},
                       'j_ubvri': {'column': [6, 7, 4, 11],
                                   'converting_method': convert_j_color}}

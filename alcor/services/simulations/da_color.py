import numpy as np


split_by_metallicities = False

ROWS_COUNT = 650
FILES_PATHS = ['cox_0524.dat',
               'cox_0570.dat',
               'cox_0593.dat',
               'cox_0609.dat',
               'cox_0632.dat',
               'cox_0659.dat',
               'cox_0705.dat',
               'cox_0767.dat',
               'cox_0837.dat',
               'cox_0877.dat']
FILES_DIR = 'da_color_tables'
MASSES = [np.array([0.524, 0.570, 0.593, 0.610, 0.632,
                    0.659, 0.705, 0.767, 0.837, 0.877])]
SEQUENCE_FILL_RULES = {'luminosity': {'column': 2},
                       'u_ubvri_absolute': {'column': 23},
                       'b_ubvri_absolute': {'column': 24},
                       'v_ubvri_absolute': {'column': 25},
                       'r_ubvri_absolute': {'column': 26},
                       'i_ubvri_absolute': {'column': 27}}

import numpy as np

ROWS_COUNT = 650
FILES_PATHS = {1: ['wd0505_z0001.trk',
                   'wd0553_z0001.trk',
                   'wd0593_z0001.trk',
                   'wd0627_z0001.trk',
                   'wd0660_z0001.trk',
                   'wd0692_z0001.trk',
                   'wd0863_z0001.trk'],
               10: ['wd0524_z001.trk',
                    'wd0570_z001.trk',
                    'wd0593_z001.trk',
                    'wd0609_z001.trk',
                    'wd0632_z001.trk',
                    'wd0659_z001.trk',
                    'wd0705_z001.trk',
                    'wd0767_z001.trk',
                    'wd0837_z001.trk',
                    'wd0877_z001.trk'],
               30: ['0524_003_sflhdiff.trk',
                    '0570_003_sflhdiff.trk',
                    '0593_003_sflhdiff.trk',
                    '0610_003_sflhdiff.trk',
                    '0632_003_sflhdiff.trk',
                    '0659_003_sflhdiff.trk',
                    '0705_003_sflhdiff.trk',
                    '1000_003_sflhdiff.trk'],
               60: ['0524_006_sflhdiff.trk',
                    '0570_006_sflhdiff.trk',
                    '0593_006_sflhdiff.trk',
                    '0610_006_sflhdiff.trk',
                    '0632_006_sflhdiff.trk',
                    '0659_006_sflhdiff.trk',
                    '0705_006_sflhdiff.trk',
                    '1000_006_sflhdiff.trk']}
FILES_FOLDER = {1:  './da_cooling_tables/Z0001',
                10: './da_cooling_tables/Z001',
                30: './da_cooling_tables/Z003',
                60: './da_cooling_tables/Z006'}
MASSES = [np.array([0.505, 0.553, 0.593, 0.627, 0.660, 0.692, 0.863]),
          np.array([0.524, 0.570, 0.593, 0.609, 0.632,
                    0.659, 0.705, 0.767, 0.837, 0.877]),
          np.array([0.524, 0.570, 0.593, 0.609,
                    0.632, 0.659, 0.705, 1.000]),
          np.array([0.524, 0.570, 0.593, 0.609,
                    0.632, 0.659, 0.705, 1.000])]
PRE_WD_LIFETIMES = [
        np.zeros(len(FILES_PATHS[1])),
        np.zeros(len(FILES_PATHS[10])),
        np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892, 0.7422, 0.4431, 0.0]),
        np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892, 0.7422, 0.4431, 0.0])]
METALLICITIES_PER_THOUSAND = [1, 10, 30, 60]
FILL_TYPES_BY_METALLICITIES = {1: 0,
                               10: 0,
                               30: 1,
                               60: 1}


def convert_effective_temperature(value: str) -> float:
    return 10. ** float(value)


def convert_cooling_time(value: str,
                         fill_type: int,
                         pre_wd_lifetime: float) -> float:
    if fill_type == 0:
        return float(value) / 1e3
    elif fill_type == 1:
        return 10. ** float(value) / 1e3 - pre_wd_lifetime

SEQUENCE_FILL_RULES = {'luminosity': {'column': 0},
                       'effective_temperature': {
                           'column': 11,
                           'converting_method': convert_effective_temperature},
                       'cooling_time': {
                           'column': [4, 8],
                           'converting_method': convert_cooling_time},
                       'surface_gravity': {'column': [11, 22]}}

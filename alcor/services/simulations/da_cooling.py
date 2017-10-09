import numpy as np

split_by_metallicities = True

ROWS_COUNT = 650
METALLICITIES_PER_THOUSAND = [1, 10, 30, 60]
FILES_PATHS = dict(zip(METALLICITIES_PER_THOUSAND,
                       [['wd0505_z0001.trk',
                         'wd0553_z0001.trk',
                         'wd0593_z0001.trk',
                         'wd0627_z0001.trk',
                         'wd0660_z0001.trk',
                         'wd0692_z0001.trk',
                         'wd0863_z0001.trk'],
                        ['wd0524_z001.trk',
                         'wd0570_z001.trk',
                         'wd0593_z001.trk',
                         'wd0609_z001.trk',
                         'wd0632_z001.trk',
                         'wd0659_z001.trk',
                         'wd0705_z001.trk',
                         'wd0767_z001.trk',
                         'wd0837_z001.trk',
                         'wd0877_z001.trk'],
                        ['0524_003_sflhdiff.trk',
                         '0570_003_sflhdiff.trk',
                         '0593_003_sflhdiff.trk',
                         '0610_003_sflhdiff.trk',
                         '0632_003_sflhdiff.trk',
                         '0659_003_sflhdiff.trk',
                         '0705_003_sflhdiff.trk',
                         '1000_003_sflhdiff.trk'],
                        ['0524_006_sflhdiff.trk',
                         '0570_006_sflhdiff.trk',
                         '0593_006_sflhdiff.trk',
                         '0610_006_sflhdiff.trk',
                         '0632_006_sflhdiff.trk',
                         '0659_006_sflhdiff.trk',
                         '0705_006_sflhdiff.trk',
                         '1000_006_sflhdiff.trk']]))
FILES_DIR = dict(zip(METALLICITIES_PER_THOUSAND,
                     ['da_cooling_tables/Z0001',
                      'da_cooling_tables/Z001',
                      'da_cooling_tables/Z003',
                      'da_cooling_tables/Z006']))
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
FILL_TYPES_BY_METALLICITIES = dict(zip(METALLICITIES_PER_THOUSAND,
                                       [0, 0, 1, 1]))


def convert_effective_temperature(value: float,
                                  fill_type: int,
                                  pre_wd_lifetime: float) -> float:
    return 10. ** value


def convert_cooling_time(value: float,
                         fill_type: int,
                         pre_wd_lifetime: float) -> float:
    if fill_type == 0:
        return value / 1e3
    elif fill_type == 1:
        return 10. ** value / 1e3 - pre_wd_lifetime

SEQUENCE_FILL_RULES = {'luminosity': {'column': [0, 0]},
                       'effective_temperature': {
                           'column': [1, 1],
                           'converting_method': convert_effective_temperature},
                       'cooling_time': {
                           'column': [4, 8],
                           'converting_method': convert_cooling_time},
                       'surface_gravity': {'column': [11, 22]}}

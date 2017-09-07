import numpy as np

ROWS_COUNT = 400
FILES_PATHS = {1: ['05047_db_Z=0.001.trk',
                   '05527_db_Z=0.001.trk',
                   '059328_db_Z=0.001.trk',
                   '062738_db_Z=0.001.trk',
                   '06602_db_Z=0.001.trk',
                   '069289_db_Z=0.001.trk',
                   '08637_db_Z=0.001.trk'],
               10: ['db_cool_0514.seq',
                    'db_cool_0530.seq',
                    'db_cool_0542.seq',
                    'db_cool_0565.seq',
                    'db_cool_0584.seq',
                    'db_cool_0609.seq',
                    'db_cool_0664.seq',
                    'db_cool_0741.seq',
                    'db_cool_0869.seq'],
               60: ['0524db_006_sflhdiff.trk',
                    '0570db_006_sflhdiff.trk',
                    '0593db_006_sflhdiff.trk',
                    '061db_006_sflhdiff.trk',
                    '0632db_006_sflhdiff.trk',
                    '0659db_006_sflhdiff.trk',
                    '070db_006_sflhdiff.trk',
                    '076db_006_sflhdiff.trk',
                    '087db_006_sflhdiff.trk']}
FILES_FOLDER = {1:  './db_cooling_tables/Z0001',
                10: './db_cooling_tables/Z001',
                60: './db_cooling_tables/Z006'}
MASSES = [np.array([0.5047, 0.5527, 0.59328, 0.62738,
                    0.6602, 0.69289, 0.8637]),
          np.array([0.514, 0.53, 0.542, 0.565, 0.584,
                    0.61, 0.664, 0.741, 0.869]),
          np.array([0.524, 0.570, 0.593, 0.61, 0.632,
                    0.659, 0.70, 0.76, 0.87])]
PRE_WD_LIFETIMES = [np.zeros(len(FILES_PATHS[1])),
                    np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892,
                              0.7422, 0.4431, 0.287, 0.114]),
                    np.array([11.117, 2.7004, 1.699, 1.2114, 0.9892,
                              0.7422, 0.4431, 0.287, 0.114])]
METALLICITIES_PER_THOUSAND = [1, 10, 60]
FILL_TYPES_BY_METALLICITIES = {1: 0,
                               10: 1,
                               60: 0}


def convert_effective_temperature(value: str,
                                  fill_type: int,
                                  pre_wd_lifetime: float) -> float:
    return 10. ** float(value)


def convert_cooling_time(value: str,
                         fill_type: int,
                         pre_wd_lifetime: float) -> float:
    if fill_type == 0:
        return 10. ** float(value) / 1e3 - pre_wd_lifetime
    elif fill_type == 1:
        return float(value)

SEQUENCE_FILL_RULES = {'luminosity': {'column': [0, 0]},
                       'effective_temperature': {
                           'column': [1, 1],
                           'converting_method': convert_effective_temperature},
                       'cooling_time': {
                           'column': [8, 2],
                           'converting_method': convert_cooling_time},
                       'surface_gravity': {'column': [22, 3]}}

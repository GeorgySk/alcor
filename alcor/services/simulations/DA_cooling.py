import csv

import numpy as np

from alcor.types import CoolingSequenceType

COOLING_SEQUENCES_FILES = {'metallicity = 0.001': ['wd0505_z0001.trk',
                                                   'wd0553_z0001.trk',
                                                   'wd0593_z0001.trk',
                                                   'wd0627_z0001.trk',
                                                   'wd0660_z0001.trk',
                                                   'wd0692_z0001.trk',
                                                   'wd0863_z0001.trk'],
                           'metallicity = 0.01': ['wd0524_z001.trk',
                                                  'wd0570_z001.trk',
                                                  'wd0593_z001.trk',
                                                  'wd0609_z001.trk',
                                                  'wd0632_z001.trk',
                                                  'wd0659_z001.trk',
                                                  'wd0705_z001.trk',
                                                  'wd0767_z001.trk',
                                                  'wd0837_z001.trk',
                                                  'wd0877_z001.trk'],
                           'metallicity = 0.03': ['0524_003_sflhdiff.trk',
                                                  '0570_003_sflhdiff.trk',
                                                  '0593_003_sflhdiff.trk',
                                                  '0610_003_sflhdiff.trk',
                                                  '0632_003_sflhdiff.trk',
                                                  '0659_003_sflhdiff.trk',
                                                  '0705_003_sflhdiff.trk',
                                                  '1000_003_sflhdiff.trk'],
                           'metallicity = 0.06': ['0524_006_sflhdiff.trk',
                                                  '0570_006_sflhdiff.trk',
                                                  '0593_006_sflhdiff.trk',
                                                  '0610_006_sflhdiff.trk',
                                                  '0632_006_sflhdiff.trk',
                                                  '0659_006_sflhdiff.trk',
                                                  '0705_006_sflhdiff.trk',
                                                  '1000_006_sflhdiff.trk']}


def initialilze_sequences() -> CoolingSequenceType:
    cooling_sequences_of_DA_WDs = {
        'metallicity = 0.001':
            {'mass': np.array(
                [0.505, 0.553, 0.593, 0.627, 0.660, 0.692, 0.863]),
                'preWD_lifetime': np.zeros(7),
                'cooling_time': np.full((7, 650), np.nan),
                'effective_temperature': np.full((7, 650), np.nan),
                'surface_gravity': np.full((7, 650), np.nan),
                'luminosity': np.full((7, 650), np.nan),
                'rows_counts': np.empty(7)},
        'metallicity = 0.01':
            {'mass': np.array(
                [0.524, 0.570, 0.593, 0.609, 0.632,
                 0.659, 0.705, 0.767, 0.837, 0.877]),
                'preWD_lifetime': np.zeros(10),
                'cooling_time': np.full((10, 650), np.nan),
                'effective_temperature': np.full((10, 650), np.nan),
                'surface_gravity': np.full((10, 650), np.nan),
                'luminosity': np.full((10, 650), np.nan),
                'rows_counts': np.empty(10)},
        'metallicity = 0.03':
            {'mass': np.array([0.524, 0.570, 0.593, 0.609,
                               0.632, 0.659, 0.705, 1.000]),
             'preWD_lifetime': np.array([11.117, 2.7004, 1.699, 1.2114,
                                         0.9892, 0.7422, 0.4431, 0.0]),
             'cooling_time': np.full((8, 650), np.nan),
             'effective_temperature': np.full((8, 650), np.nan),
             'surface_gravity': np.full((8, 650), np.nan),
             'luminosity': np.full((8, 650), np.nan),
             'rows_counts': np.empty(8)
             },
        'metallicity = 0.06':
            {'mass': np.array([0.524, 0.570, 0.593, 0.609,
                               0.632, 0.659, 0.705, 1.000]),
             'preWD_lifetime': np.array([11.117, 2.7004, 1.699, 1.2114,
                                         0.9892, 0.7422, 0.4431, 0.0]),
             'cooling_time': np.full((8, 650), np.nan),
             'effective_temperature': np.full((8, 650), np.nan),
             'surface_gravity': np.full((8, 650), np.nan),
             'luminosity': np.full((8, 650), np.nan),
             'rows_counts': np.empty(8)
             }
    }

    for filename_index, filename in COOLING_SEQUENCES_FILES[
            'metallicity = 0.001']:
        cooling_sequences = cooling_sequences_of_DA_WDs['metallicity = 0.001']
        cooling_time = cooling_sequences['cooling_time']
        effective_temperature = cooling_sequences['effective_temperature']
        surface_gravity = cooling_sequences['surface_gravity']
        luminosity = cooling_sequences['luminosity']
        rows_counts = cooling_sequences['rows_counts']

        with open(filename, 'r') as file:
            filereader = csv.reader(file,
                                    delimiter=' ',
                                    skipinitialspace=True)

            rows_counts[filename_index] = sum(1 for row in filereader)

            for row_index, row in enumerate(filereader):
                cooling_time[filename_index, row_index] = (float(row[5])
                                                           / 1000.0)
                effective_temperature[filename_index, row_index] = (
                    10. ** float(row[1]))
                surface_gravity[filename_index, row_index] = float(row[11])
                luminosity[filename_index, row_index] = float(row[0])

    for filename_index, filename in COOLING_SEQUENCES_FILES[
            'metallicity = 0.01']:
        cooling_sequences = cooling_sequences_of_DA_WDs['metallicity = 0.01']
        cooling_time = cooling_sequences['cooling_time']
        effective_temperature = cooling_sequences['effective_temperature']
        surface_gravity = cooling_sequences['surface_gravity']
        luminosity = cooling_sequences['luminosity']
        rows_counts = cooling_sequences['rows_counts']

        with open(filename, 'r') as file:
            filereader = csv.reader(file,
                                    delimiter=' ',
                                    skipinitialspace=True)

            rows_counts[filename_index] = sum(1 for row in filereader)

            for row_index, row in enumerate(filereader):
                cooling_time[filename_index, row_index] = (float(row[5])
                                                           / 1000.0)
                effective_temperature[filename_index, row_index] = (
                    10. ** float(row[1]))
                surface_gravity[filename_index, row_index] = float(row[11])
                luminosity[filename_index, row_index] = float(row[0])

    for filename_index, filename in COOLING_SEQUENCES_FILES[
            'metallicity = 0.03']:
        cooling_sequences = cooling_sequences_of_DA_WDs['metallicity = 0.03']
        preWD_lifetime = cooling_sequences['preWD_lifetime']
        cooling_time = cooling_sequences['cooling_time']
        effective_temperature = cooling_sequences['effective_temperature']
        surface_gravity = cooling_sequences['surface_gravity']
        luminosity = cooling_sequences['luminosity']
        rows_counts = cooling_sequences['rows_counts']

        with open(filename, 'r') as file:
            filereader = csv.reader(file,
                                    delimiter=' ',
                                    skipinitialspace=True)

            rows_counts[filename_index] = sum(1 for row in filereader)

            for row_index, row in enumerate(filereader):
                cooling_time[filename_index, row_index] = (
                    10. ** float(row[8]) / 1000.
                    - preWD_lifetime[filename_index])
                effective_temperature[filename_index, row_index] = (
                    10. ** float(row[1]))
                surface_gravity[filename_index, row_index] = float(row[22])
                luminosity[filename_index, row_index] = float(row[0])

    for filename_index, filename in COOLING_SEQUENCES_FILES[
            'metallicity = 0.06']:
        cooling_sequences = cooling_sequences_of_DA_WDs['metallicity = 0.06']
        preWD_lifetime = cooling_sequences['preWD_lifetime']
        cooling_time = cooling_sequences['cooling_time']
        effective_temperature = cooling_sequences['effective_temperature']
        surface_gravity = cooling_sequences['surface_gravity']
        luminosity = cooling_sequences['luminosity']
        rows_counts = cooling_sequences['rows_counts']

        with open(filename, 'r') as file:
            filereader = csv.reader(file,
                                    delimiter=' ',
                                    skipinitialspace=True)

            rows_counts[filename_index] = sum(1 for row in filereader)

            for row_index, row in enumerate(filereader):
                cooling_time[filename_index, row_index] = (
                    10. ** float(row[8]) / 1000.
                    - preWD_lifetime[filename_index])
                effective_temperature[filename_index, row_index] = (
                    10. ** float(row[1]))
                surface_gravity[filename_index, row_index] = float(row[22])
                luminosity[filename_index, row_index] = float(row[0])

    return cooling_sequences_of_DA_WDs

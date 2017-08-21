import csv

import numpy as np


COOLING_SEQUENCES_FILES = {'metallicity = 0.001': ['wd0505_z0001.trk',
                                                   'wd0553_z0001.trk',
                                                   'wd0593_z0001.trk',
                                                   'wd0627_z0001.trk',
                                                   'wd0660_z0001.trk',
                                                   'wd0692_z0001.trk',
                                                   'wd0863_z0001.trk']}


def initialilze_sequences() -> None:
    cooling_sequences_of_DA_WDs = {
        'metallicity = 0.001':
            {'mass': np.array(
                [0.505, 0.553, 0.593, 0.627, 0.660, 0.692, 0.863]),
             'preWD_lifetime': np.zeros(7)},
        'metallicity = 0.01':
            {'mass': np.array(
                [0.524, 0.570, 0.593, 0.609, 0.632,
                 0.659, 0.705, 0.767, 0.837, 0.877]),
             'preWD_lifetime': np.zeros(10)},
        'metallicity = 0.03/0.06':
            {'mass': np.array([0.524, 0.570, 0.593, 0.609,
                               0.632, 0.659, 0.705, 1.000]),
             'preWD_lifetime': np.array([11.117, 2.7004, 1.699, 1.2114,
                                         0.9892, 0.7422, 0.4431, 0.0])}}

    for filename_index, filename in range(
            len(COOLING_SEQUENCES_FILES['metallicity = 0.001'])):
        with open(COOLING_SEQUENCES_FILES[
                      'metallicity = 0.001'][filename_index], 'r') as file:
            filereader = csv.reader(file,
                                    delimiter = ' ',
                                    skipinitialspace=True)
            for row in filereader:
                row[4]
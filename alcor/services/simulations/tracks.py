import os
import posixpath
from collections import OrderedDict
from contextlib import (contextmanager,
                        closing)
from functools import partial
from typing import (Dict,
                    Tuple,
                    List)

import h5py
import pandas as pd


join_group = posixpath.join

COLORS = ['color_u',
          'color_b',
          'color_v',
          'color_r',
          'color_i',
          'color_j']


def fill_cooling_tracks(cooling_tracks: Dict[int, Dict],
                        *,
                        file: h5py.File) -> None:
    for metallicity in cooling_tracks.keys():
        cooling_tracks_by_metallicity = cooling_tracks[metallicity]

        metallicity_group = str(metallicity)
        masses_strings = sort_mass_indexes(file[metallicity_group])

        for mass in masses_strings:
            mass_group = join_group(metallicity_group, mass)
            tracks = OrderedDict((key, file[join_group(mass_group,
                                                       key)])
                                 for key in ['cooling_time',
                                             'effective_temperature',
                                             'luminosity'])
            cooling_tracks_by_metallicity[int(mass)] = pd.DataFrame(tracks)


def read_cooling(path: str,
                 metallicities: Tuple[int, ...]
                 ) -> Dict[int, Dict[int, pd.DataFrame]]:
    with open_hdf5(path) as file:
        cooling_tracks_by_metallicities = OrderedDict(
                (metallicity, OrderedDict())
                for metallicity in metallicities)
        fill_cooling_tracks(cooling_tracks_by_metallicities,
                            file=file)
        return cooling_tracks_by_metallicities


# TODO: all the paths must be passed from upper level module
read_da_cooling = partial(read_cooling,
                          path='input_data/da_cooling.hdf5',
                          metallicities=(1, 10, 30, 60))
read_db_cooling = partial(read_cooling,
                          path='input_data/db_cooling.hdf5',
                          metallicities=(1, 10, 60))


def fill_table(*,
               file: h5py.File,
               interest_parameters: List[str]) -> Dict[int, pd.DataFrame]:
    for mass_group in file:
        tracks = OrderedDict((parameter, file[join_group(mass_group,
                                                         parameter)])
                             for parameter in interest_parameters)
        yield int(mass_group), pd.DataFrame(tracks)


def read_table(path: str,
               *,
               interest_parameters: List[str]) -> Dict[int, pd.DataFrame]:
    with open_hdf5(path) as file:
        return OrderedDict(fill_table(
                file=file,
                interest_parameters=interest_parameters))


read_da_db_colors = partial(read_table,
                            interest_parameters=['luminosity'] + COLORS)

read_da_colors = partial(read_da_db_colors,
                         path='input_data/da_colors.hdf5')

read_db_colors = partial(read_da_db_colors,
                         path='input_data/db_colors.hdf5')

read_one_tables = partial(read_table,
                          path='input_data/one_wds_tracks.hdf5',
                          interest_parameters=(['luminosity',
                                                'cooling_time',
                                                'effective_temperature']
                                               + COLORS))


@contextmanager
def open_hdf5(path: str) -> h5py.File:
    base_dir = os.path.dirname(__file__)
    path = os.path.join(base_dir, path)
    with closing(h5py.File(path,
                           mode='r')) as file:
        yield file


def sort_mass_indexes(indexes: List[str]) -> List[str]:
    integer_indexes = sorted(map(int, indexes))
    return list(map(str, integer_indexes))

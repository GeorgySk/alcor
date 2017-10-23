import os
import posixpath
from contextlib import (contextmanager,
                        closing)
from functools import partial
from typing import (Dict,
                    Tuple,
                    List)

import h5py
import pandas as pd


join_group = posixpath.join


def read_cooling(path: str,
                 metallicities: Tuple[int]
                 ) -> Dict[int, Dict[int, pd.DataFrame]]:
    with open_hdf5(path) as file:
        cooling_tracks_by_metallicities = {metallicity: {}
                                           for metallicity in metallicities}
        fill_cooling_tracks(cooling_tracks_by_metallicities,
                            file=file)
        return cooling_tracks_by_metallicities


read_da_cooling = partial(read_cooling,
                          path='input_data/da_cooling.hdf5',
                          metallicities=(1, 10, 30, 60))
read_db_cooling = partial(read_cooling,
                          path='input_data/db_cooling.hdf5',
                          metallicities=(1, 10, 60))


def read_colors(path: str) -> Dict[int, pd.DataFrame]:
    with open_hdf5(path) as file:
        return fill_colors(file=file)


read_da_colors = partial(read_colors,
                         path='input_data/da_colors.hdf5')

read_db_colors = partial(read_colors,
                         path='input_data/db_colors.hdf5')


def read_one_tables(path: str = 'input_data/one_wds_tracks.hdf5'
                    ) -> Dict[int, pd.DataFrame]:
    with open_hdf5(path) as file:
        return fill_one_table(file=file)


def fill_cooling_tracks(cooling_tracks: Dict[int, Dict],
                        *,
                        file: h5py.File) -> None:
    for metallicity in cooling_tracks.keys():
        metallicity_group = str(metallicity)
        masses = sort_mass_indexes(indexes=file[metallicity_group])

        for mass in masses:
            mass_group = join_group(metallicity_group, mass)
            cooling_tracks_by_metallicity = cooling_tracks[metallicity]
            cooling_tracks_by_metallicity[int(mass)] = pd.DataFrame(
                    dict(cooling_time=file[join_group(mass_group,
                                                      'cooling_time')],
                         effective_temperature=file[join_group(
                                 mass_group, 'effective_temperature')],
                         luminosity=file[join_group(mass_group,
                                                    'luminosity')]))


def fill_colors(*,
                file: h5py.File) -> Dict[int, pd.DataFrame]:
    color_table = {}
    for mass_group in file:
        color_table[int(mass_group)] = pd.DataFrame(dict(
                luminosity=file[join_group(mass_group, 'luminosity')],
                color_u=file[join_group(mass_group, 'color_u')],
                color_b=file[join_group(mass_group, 'color_b')],
                color_v=file[join_group(mass_group, 'color_v')],
                color_r=file[join_group(mass_group, 'color_r')],
                color_i=file[join_group(mass_group, 'color_i')],
                color_j=file[join_group(mass_group, 'color_j')]))
    return color_table


def fill_one_table(*,
                   file: h5py.File) -> Dict[int, pd.DataFrame]:
    table = {}
    for mass_group in file:
        table[int(mass_group)] = pd.DataFrame(dict(
                luminosity=file[join_group(mass_group, 'luminosity')],
                cooling_time=file[join_group(mass_group, 'cooling_time')],
                effective_temperature=file[join_group(
                        mass_group, 'effective_temperature')],
                color_u=file[join_group(mass_group, 'color_u')],
                color_b=file[join_group(mass_group, 'color_b')],
                color_v=file[join_group(mass_group, 'color_v')],
                color_r=file[join_group(mass_group, 'color_r')],
                color_i=file[join_group(mass_group, 'color_i')],
                color_j=file[join_group(mass_group, 'color_j')]))
    return table


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

import os
import posixpath
from collections import OrderedDict
from contextlib import (contextmanager,
                        closing)
from typing import (Dict,
                    List,
                    Tuple)

import h5py
import pandas as pd


join_group = posixpath.join


def fill_cooling_tracks(cooling_tracks: Dict[int, Dict],
                        *,
                        file: h5py.File,
                        interest_parameters: Tuple[str, ...]) -> None:
    for metallicity in cooling_tracks.keys():
        cooling_tracks_by_metallicity = cooling_tracks[metallicity]

        metallicity_group = str(metallicity)
        masses_strings = sort_mass_indexes(file[metallicity_group])

        for mass_string in masses_strings:
            mass_group = join_group(metallicity_group, mass_string)
            tracks = OrderedDict((parameter, file[join_group(mass_group,
                                                             parameter)])
                                 for parameter in interest_parameters)
            cooling_tracks_by_metallicity[int(mass_string)] = (
                pd.DataFrame(tracks))


def read_cooling(path: str,
                 metallicities: Tuple[int, ...],
                 interest_parameters: Tuple[str, ...]
                 ) -> Dict[int, Dict[int, pd.DataFrame]]:
    with open_hdf5(path) as file:
        cooling_tracks_by_metallicities = OrderedDict(
                (metallicity, OrderedDict())
                for metallicity in metallicities)
        fill_cooling_tracks(cooling_tracks_by_metallicities,
                            file=file,
                            interest_parameters=interest_parameters)
        return cooling_tracks_by_metallicities


def fill_table(*,
               file: h5py.File,
               interest_parameters: Tuple[str, ...]
               ) -> Dict[int, pd.DataFrame]:
    for mass_group in file:
        tracks = {parameter: file[join_group(mass_group, parameter)]
                  for parameter in interest_parameters}
        yield int(mass_group), pd.DataFrame(tracks)


def read_table(path: str,
               *,
               interest_parameters: Tuple[str, ...]
               ) -> Dict[int, pd.DataFrame]:
    with open_hdf5(path) as file:
        return OrderedDict(fill_table(
                file=file,
                interest_parameters=interest_parameters))


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

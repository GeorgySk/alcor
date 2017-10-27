import os
from functools import partial
from typing import List

import pytest

from alcor.services.simulations import tracks
from tests.utils import (fort_files_lengths,
                         fort_files_by_metallicities_lengths,
                         tracks_lengths,
                         tracks_by_metallicities_lengths)

FORT_FILES_PATH = os.path.abspath('tests/tables')


@pytest.fixture(scope='function')
def da_cooling_tracks_lengths() -> List[int]:
    cooling_tracks_by_metallicities = tracks.read_da_cooling()
    return list(tracks_by_metallicities_lengths(
            cooling_tracks=cooling_tracks_by_metallicities))


@pytest.fixture(scope='function')
def db_cooling_tracks_lengths() -> List[int]:
    cooling_tracks_by_metallicities = tracks.read_db_cooling()
    return list(tracks_by_metallicities_lengths(
            cooling_tracks=cooling_tracks_by_metallicities))


@pytest.fixture(scope='function')
def da_colors_tracks_lengths() -> List[int]:
    characteristics = tracks.read_da_colors()
    return list(tracks_lengths(characteristics))


@pytest.fixture(scope='function')
def db_colors_tracks_lengths() -> List[int]:
    characteristics = tracks.read_db_colors()
    return list(tracks_lengths(characteristics))


@pytest.fixture(scope='function')
def one_tables_lengths() -> List[int]:
    characteristics = tracks.read_one_tables()
    return list(tracks_lengths(characteristics))


@pytest.fixture(scope='function')
def da_cooling_fort_files_lengths() -> List[int]:
    da_cooling_fort_links = {1: range(11, 18),
                             10: range(21, 31),
                             30: range(31, 39),
                             60: range(41, 49)}
    return list(fort_files_by_metallicities_lengths(
            fort_links=da_cooling_fort_links,
            base_dir=FORT_FILES_PATH))


@pytest.fixture(scope='function')
def db_cooling_fort_files_lengths() -> List[int]:
    db_cooling_fort_links = {1: range(91, 98),
                             10: range(101, 110),
                             60: range(111, 120)}
    return list(fort_files_by_metallicities_lengths(
            fort_links=db_cooling_fort_links,
            base_dir=FORT_FILES_PATH))


forts_lengths = partial(fort_files_lengths,
                        base_dir=FORT_FILES_PATH)


@pytest.fixture(scope='function')
def da_colors_fort_files_lengths() -> List[int]:
    da_colors_fort_links = range(61, 71)
    return list(forts_lengths(fort_links=da_colors_fort_links))


@pytest.fixture(scope='function')
def db_colors_fort_files_lengths() -> List[int]:
    db_colors_fort_links = range(132, 139)
    return list(forts_lengths(fort_links=db_colors_fort_links))


@pytest.fixture(scope='function')
def one_tables_fort_files_lengths() -> List[int]:
    one_tables_fort_links = range(121, 127)
    return list(forts_lengths(fort_links=one_tables_fort_links))

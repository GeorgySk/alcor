import os
from typing import (Dict,
                    List)

import pandas as pd
import pytest

from alcor.services.simulations import tracks
from tests.utils import (fort_files_by_metallicities_lengths,
                         fort_files_lengths)

FORT_FILES_PATH = os.path.abspath('/alcor/tests/tables')


@pytest.fixture
def da_cooling_tracks_lengths() -> List[int]:
    cooling_tracks_by_metallicities = tracks.read_da_cooling()
    return tracks_by_metallicities_lengths(
            cooling_tracks=cooling_tracks_by_metallicities)


@pytest.fixture
def db_cooling_tracks_lengths() -> List[int]:
    cooling_tracks_by_metallicities = tracks.read_db_cooling()
    return tracks_by_metallicities_lengths(
            cooling_tracks=cooling_tracks_by_metallicities)


@pytest.fixture
def da_colors_tracks_lengths() -> List[int]:
    characteristics = tracks.read_da_colors()
    return tracks_lengths(sequences=characteristics)


@pytest.fixture
def db_colors_tracks_lengths() -> List[int]:
    characteristics = tracks.read_db_colors()
    return tracks_lengths(sequences=characteristics)


@pytest.fixture
def one_tables_lengths() -> List[int]:
    characteristics = tracks.read_one_tables()
    return tracks_lengths(sequences=characteristics)


@pytest.fixture
def da_cooling_fort_files_lengths() -> List[int]:
    da_cooling_fort_links = {1: range(11, 18),
                             10: range(21, 31),
                             30: range(31, 39),
                             60: range(41, 49)}
    return list(fort_files_by_metallicities_lengths(
            fort_links=da_cooling_fort_links,
            base_dir=FORT_FILES_PATH))


@pytest.fixture
def db_cooling_fort_files_lengths() -> List[int]:
    db_cooling_fort_links = {1: range(91, 98),
                             10: range(101, 110),
                             60: range(111, 120)}
    return list(fort_files_by_metallicities_lengths(
            fort_links=db_cooling_fort_links,
            base_dir=FORT_FILES_PATH))


@pytest.fixture
def da_colors_fort_files_lengths() -> List[int]:
    da_colors_fort_links = range(61, 71)
    return list(fort_files_lengths(
            fort_links=da_colors_fort_links,
            base_dir=FORT_FILES_PATH))


@pytest.fixture
def db_colors_fort_files_lengths() -> List[int]:
    db_colors_fort_links = range(132, 139)
    return list(fort_files_lengths(
            fort_links=db_colors_fort_links,
            base_dir=FORT_FILES_PATH))


@pytest.fixture
def one_fort_files_lengths() -> List[int]:
    one_tables_fort_links = range(121, 127)
    return list(fort_files_lengths(
            fort_links=one_tables_fort_links,
            base_dir=FORT_FILES_PATH))


@pytest.fixture
def tracks_by_metallicities_lengths(
        cooling_tracks: Dict[int, Dict[int, pd.DataFrame]]) -> List[int]:
    lengths = []
    for metallicity, sequences in cooling_tracks.items():
        lengths.extend(tracks_lengths(sequences=sequences))

    return lengths


@pytest.fixture
def tracks_lengths(sequences: Dict[int, pd.DataFrame]) -> List[int]:
    lengths = []
    for mass, sequence in sorted(sequences.items()):
        lengths.append(sequence.shape[0])

    return lengths

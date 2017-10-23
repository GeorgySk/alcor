import os
from typing import (Iterator,
                    Dict,
                    List)

import pandas as pd
import pytest

from alcor.services.simulations import tracks


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
def da_cooling_fort_files_lengths():
    da_cooling_fort_links = {1: range(11, 18),
                             10: range(21, 31),
                             30: range(31, 39),
                             60: range(41, 49)}
    return fort_files_by_metallicities_lengths(
            fort_links=da_cooling_fort_links,
            base_dir=FORT_FILES_PATH)


@pytest.fixture
def db_cooling_fort_files_lengths():
    db_cooling_fort_links = {1: range(91, 98),
                             10: range(101, 110),
                             60: range(111, 120)}
    return fort_files_by_metallicities_lengths(
            fort_links=db_cooling_fort_links,
            base_dir=FORT_FILES_PATH)


@pytest.fixture
def da_colors_fort_files_lengths():
    da_colors_fort_links = range(61, 71)
    return fort_files_lengths(
            fort_links=da_colors_fort_links,
            base_dir=FORT_FILES_PATH)


@pytest.fixture
def db_colors_fort_files_lengths():
    db_colors_fort_links = range(132, 139)
    return fort_files_lengths(
            fort_links=db_colors_fort_links,
            base_dir=FORT_FILES_PATH)


@pytest.fixture
def one_fort_files_lengths():
    one_tables_fort_links = range(121, 127)
    return fort_files_lengths(
            fort_links=one_tables_fort_links,
            base_dir=FORT_FILES_PATH)


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
    for mass, sequence in sequences.items():
        lengths.append(sequence.shape[0])

    return lengths


@pytest.fixture
def fort_files_by_metallicities_lengths(fort_links: Dict[int, Iterator],
                                        *,
                                        base_dir: str) -> List[int]:
    lengths = []
    for metallicity, fort_links_range in fort_links.items():
        lengths.extend(fort_files_lengths(fort_links=fort_links_range,
                                          base_dir=base_dir))
    return lengths


@pytest.fixture
def fort_files_lengths(fort_links: Iterator,
                       *,
                       base_dir: str) -> List[int]:
    lengths = []
    fort_link_dirs = [os.path.join(base_dir,
                                   'fort_files/fort.' + str(fort_link))
                      for fort_link in fort_links]
    for fort_link_dir in fort_link_dirs:
        lengths.append(sum(1 for line in open(fort_link_dir)))

    return lengths

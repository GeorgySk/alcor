import os
from typing import (Iterator,
                    Dict,
                    List)

import pandas as pd

from alcor.services.simulations import tracks


def test_tracks() -> None:
    da_cooling_tracks_by_metallicities = tracks.read_da_cooling()
    db_cooling_tracks_by_metallicities = tracks.read_db_cooling()
    da_colors = tracks.read_da_colors()
    db_colors = tracks.read_db_colors()
    one_tables = tracks.read_one_tables()

    da_cooling_fort_links = {1: range(11, 18),
                             10: range(21, 31),
                             30: range(31, 39),
                             60: range(41, 49)}
    db_cooling_fort_links = {1: range(91, 98),
                             10: range(101, 110),
                             60: range(111, 120)}
    da_colors_fort_links = range(61, 71)
    db_colors_fort_links = range(132, 139)
    one_tables_fort_links = range(121, 127)

    base_dir = os.path.dirname(__file__)

    da_cooling_tracks_lengths = tracks_by_metallicities_lengths(
            cooling_tracks=da_cooling_tracks_by_metallicities)
    db_cooling_tracks_lengths = tracks_by_metallicities_lengths(
            cooling_tracks=db_cooling_tracks_by_metallicities)

    da_cooling_fort_files_lengths = fort_files_by_metallicities_lengths(
            fort_links=da_cooling_fort_links,
            base_dir=base_dir)
    db_cooling_fort_files_lengths = fort_files_by_metallicities_lengths(
            fort_links=db_cooling_fort_links,
            base_dir=base_dir)

    da_colors_tracks_lengths = tracks_lengths(sequences=da_colors)
    db_colors_tracks_lengths = tracks_lengths(sequences=db_colors)

    da_colors_fort_files_lengths = fort_files_lengths(
            fort_links=da_colors_fort_links,
            base_dir=base_dir)
    db_colors_fort_files_lengths = fort_files_lengths(
            fort_links=db_colors_fort_links,
            base_dir=base_dir)

    one_tracks_lengths = tracks_lengths(sequences=one_tables)
    one_fort_files_lengths = fort_files_lengths(
            fort_links=one_tables_fort_links,
            base_dir=base_dir)

    assert da_cooling_tracks_lengths == da_cooling_fort_files_lengths
    assert db_cooling_tracks_lengths == db_cooling_fort_files_lengths
    assert da_colors_tracks_lengths == da_colors_fort_files_lengths
    assert db_colors_tracks_lengths == db_colors_fort_files_lengths
    assert one_tracks_lengths == one_fort_files_lengths


def tracks_by_metallicities_lengths(
        cooling_tracks: Dict[int, Dict[int, pd.DataFrame]]) -> List[int]:
    lengths = []
    for metallicity, sequences in cooling_tracks.items():
        lengths.extend(tracks_lengths(sequences=sequences))

    return lengths


def tracks_lengths(sequences: Dict[int, pd.DataFrame]) -> List[int]:
    lengths = []
    for mass, sequence in sequences.items():
        lengths.append(sequence.shape[0])

    return lengths


def fort_files_by_metallicities_lengths(fort_links: Dict[int, Iterator],
                                        *,
                                        base_dir: str) -> List[int]:
    lengths = []
    for metallicity, fort_links_range in fort_links.items():
        lengths.extend(fort_files_lengths(fort_links=fort_links_range,
                                          base_dir=base_dir))
    return lengths


def fort_files_lengths(fort_links: Iterator,
                       *,
                       base_dir: str) -> List[int]:
    lengths = []
    for fort_link in fort_links:
        fort_link_dir = os.path.join(base_dir,
                                     './fort_files/fort.' + str(fort_link))
        lengths.append(sum(1 for line in open(fort_link_dir)))

    return lengths


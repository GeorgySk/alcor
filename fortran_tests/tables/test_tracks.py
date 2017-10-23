from typing import List


def test_tracks(da_cooling_tracks_lengths: List[int],
                db_cooling_tracks_lengths: List[int],
                da_colors_tracks_lengths: List[int],
                db_colors_tracks_lengths: List[int],
                one_tables_lengths: List[int],
                da_cooling_fort_files_lengths: List[int],
                db_cooling_fort_files_lengths: List[int],
                da_colors_fort_files_lengths: List[int],
                db_colors_fort_files_lengths: List[int],
                one_fort_files_lengths: List[int]) -> None:
    assert da_cooling_tracks_lengths == da_cooling_fort_files_lengths
    assert db_cooling_tracks_lengths == db_cooling_fort_files_lengths
    assert da_colors_tracks_lengths == da_colors_fort_files_lengths
    assert db_colors_tracks_lengths == db_colors_fort_files_lengths
    assert one_tables_lengths == one_fort_files_lengths


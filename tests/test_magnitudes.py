from typing import (Dict,
                    Tuple,
                    List)

import math
import numpy as np
import pandas as pd

from alcor.services.simulations.magnitudes import (estimate_at,
                                                   estimated_interest_value,
                                                   calculate_index,
                                                   extrapolating_by_grid,
                                                   get_min_metallicity_index,
                                                   generate_spectral_types,
                                                   extrapolate_interest_value,
                                                   interpolate_interest_value,
                                                   estimate_color,
                                                   estimate_by_mass,
                                                   estimate_by_metallicities,
                                                   assign_estimated_values)


def test_estimate_at(float_value: float,
                     floats_tuple: Tuple[float, float],
                     other_floats_tuple: Tuple[float, float]) -> None:
    estimated_value = estimate_at(float_value,
                                  x=floats_tuple,
                                  y=other_floats_tuple)
    other_estimated_value = estimate_at(float_value,
                                        x=floats_tuple,
                                        y=(float_value, float_value))

    assert isinstance(estimated_value, float)
    assert math.isclose(float_value, other_estimated_value)


def test_estimated_interest_value(cooling_time: float,
                                  cooling_time_grid: np.ndarray,
                                  interest_parameter_grid: np.ndarray,
                                  row_index: int) -> None:
    value = estimated_interest_value(
            cooling_time=cooling_time,
            cooling_time_grid=cooling_time_grid,
            interest_parameter_grid=interest_parameter_grid,
            row_index=row_index)

    assert isinstance(value, float)


def test_calculate_index(float_value: float,
                         grid: np.ndarray) -> None:
    index = calculate_index(float_value,
                            grid=grid)

    assert isinstance(index, int)
    assert (index == -2) or (index >= 0)


def test_extrapolating_by_grid(cooling_time: float,
                               cooling_time_grid: np.ndarray) -> None:
    result = extrapolating_by_grid(cooling_time,
                                   cooling_time_grid=cooling_time_grid)

    assert isinstance(result, bool)


def test_get_min_metallicity_index(metallicity: float,
                                   grid_metallicities: List[float]) -> None:
    metallicity_index = get_min_metallicity_index(
            metallicity,
            grid_metallicities=grid_metallicities)

    assert isinstance(metallicity_index, int)


def test_generate_spectral_types(db_to_da_fraction: float,
                                 size: int) -> None:
    spectral_types = generate_spectral_types(
            db_to_da_fraction=db_to_da_fraction,
            size=size)

    assert isinstance(spectral_types, np.ndarray)


def test_extrapolate_interest_value(
        mass: float,
        cooling_time: float,
        greater_mass_cooling_time_grid: np.ndarray,
        greater_mass_interest_parameter_grid: np.ndarray,
        lesser_mass_cooling_time_grid: np.ndarray,
        lesser_mass_interest_parameter_grid: np.ndarray,
        min_mass: float,
        max_mass: float,
        min_row_index: int,
        max_row_index: int) -> None:
    value = extrapolate_interest_value(
            mass=mass,
            cooling_time=cooling_time,
            greater_mass_cooling_time_grid=greater_mass_cooling_time_grid,
            greater_mass_interest_parameter_grid=(
                greater_mass_interest_parameter_grid),
            lesser_mass_cooling_time_grid=lesser_mass_cooling_time_grid,
            lesser_mass_interest_parameter_grid=(
                lesser_mass_interest_parameter_grid),
            min_mass=min_mass,
            max_mass=max_mass,
            min_row_index=min_row_index,
            max_row_index=max_row_index)

    assert isinstance(value, float)


def test_interpolate_interest_value(
        mass: float,
        cooling_time: float,
        greater_mass_cooling_time_grid: np.ndarray,
        greater_mass_interest_parameter_grid: np.ndarray,
        lesser_mass_cooling_time_grid: np.ndarray,
        lesser_mass_interest_parameter_grid: np.ndarray,
        min_mass: float,
        max_mass: float,
        min_row_index: int,
        max_row_index: int) -> None:
    value = interpolate_interest_value(
            mass=mass,
            cooling_time=cooling_time,
            greater_mass_cooling_time_grid=greater_mass_cooling_time_grid,
            greater_mass_interest_parameter_grid=(
                greater_mass_interest_parameter_grid),
            lesser_mass_cooling_time_grid=lesser_mass_cooling_time_grid,
            lesser_mass_interest_parameter_grid=(
                lesser_mass_interest_parameter_grid),
            min_mass=min_mass,
            max_mass=max_mass,
            min_row_index=min_row_index,
            max_row_index=max_row_index)

    assert isinstance(value, float)


def test_estimate_color(star_series: pd.Series,
                        color_table: Dict[int, pd.DataFrame],
                        color: str) -> None:
    color = estimate_color(star_series,
                           color_table=color_table,
                           color=color)

    isinstance(color, float)


def test_estimate_by_mass(star_series: pd.Series,
                          tracks: Dict[int, pd.DataFrame],
                          interest_parameter: str) -> None:
    parameter = estimate_by_mass(star_series,
                                 tracks=tracks,
                                 interest_parameter=interest_parameter)

    isinstance(parameter, float)


def test_estimate_by_metallicities(
        star_series: pd.Series,
        metallicity_grid: List[float],
        cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        interest_parameter: str) -> None:
    parameter = estimate_by_metallicities(
            star_series,
            metallicity_grid=metallicity_grid,
            cooling_sequences=cooling_sequences,
            interest_parameter=interest_parameter)

    isinstance(parameter, float)


def test_assign_estimated_values(
        stars_without_luminosity: pd.DataFrame,
        da_cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        da_color_table: Dict[int, pd.DataFrame],
        db_cooling_sequences: Dict[int, Dict[int, pd.DataFrame]],
        db_color_table: Dict[int, pd.DataFrame],
        one_color_table: Dict[int, pd.DataFrame]) -> None:
    stars_without_luminosity = assign_estimated_values(
            stars_without_luminosity,
            da_cooling_sequences=da_cooling_sequences,
            da_color_table=da_color_table,
            db_cooling_sequences=db_cooling_sequences,
            db_color_table=db_color_table,
            one_color_table=one_color_table)

    isinstance(stars_without_luminosity, pd.DataFrame)

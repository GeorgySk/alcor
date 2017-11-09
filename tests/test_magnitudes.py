from typing import (Tuple,
                    List)

import numpy as np

from alcor.services.simulations.magnitudes import (estimate_at,
                                                   estimated_interest_value,
                                                   calculate_index,
                                                   extrapolating_by_grid,
                                                   get_min_metallicity_index,
                                                   generate_spectral_types)


def test_estimate_at(float_value: float,
                     floats_tuple: Tuple[float, float],
                     other_floats_tuple: Tuple[float, float]) -> None:
    estimated_value = estimate_at(float_value,
                                  x=floats_tuple,
                                  y=other_floats_tuple)

    assert isinstance(estimated_value, float)


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
    res = extrapolating_by_grid(cooling_time,
                                cooling_time_grid=cooling_time_grid)

    assert isinstance(res, bool)


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

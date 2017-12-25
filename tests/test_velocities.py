import numpy as np

from alcor.services.simulations.velocities import rotate_vectors


def test_rotate_vectors(*,
                        x_values: np.ndarray,
                        y_values: np.ndarray,
                        angles: np.ndarray) -> None:
    rotated_x_values, rotated_y_values = rotate_vectors(x_values=x_values,
                                                        y_values=y_values,
                                                        angles=angles)

    assert isinstance(rotated_x_values, np.ndarray)
    assert isinstance(rotated_y_values, np.ndarray)
    assert rotated_x_values.size == x_values.size
    assert rotated_y_values.size == y_values.size

import numpy as np
import pandas as pd

PC_PER_KPC = 1e3


def assign_proper_motions(stars: pd.DataFrame,
                          *,
                          kappa: float = 4.74,
                          sin_longitude: np.ndarray,
                          cos_longitude: np.ndarray,
                          sin_latitude: np.ndarray,
                          cos_latitude: np.ndarray) -> None:
    """
    :param stars:
    :param kappa: km/sec in one au/year
    :param sin_longitude:
    :param cos_longitude:
    :param sin_latitude:
    :param cos_latitude:
    :return:
    """
    velocities_by_proper_motion = 1. / (kappa * stars['distance'] * PC_PER_KPC)

    stars['proper_motion_in_longitude'] = (
        velocities_by_proper_motion * (- stars['u_velocity'] * sin_longitude
                                       + stars['v_velocity'] * cos_longitude))
    stars['proper_motion_in_latitude'] = velocities_by_proper_motion * (
        - stars['u_velocity'] * cos_longitude * sin_latitude
        - stars['v_velocity'] * sin_latitude * sin_longitude
        + stars['w_velocity'] * cos_latitude)
    stars['radial_velocity'] = (
        stars['u_velocity'] * cos_latitude * cos_longitude
        + stars['v_velocity'] * cos_latitude * sin_latitude
        + stars['w_velocity'] * sin_latitude)

    stars['proper_motion'] = np.sqrt(stars['proper_motion_in_longitude'] ** 2
                                     + stars['proper_motion_in_latitude'] ** 2)

from decimal import Decimal

import sys
from hypothesis import strategies

from alcor.models import Star

stars = strategies.builds(
    Star,
    # group_id
    strategies.text(min_size=6),
    # luminosity
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # proper_motion
    strategies.decimals(min_value=Decimal(sys.float_info.epsilon),
                        allow_nan=False,
                        allow_infinity=False),
    # proper_motion_component_b
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # proper_motion_component_l
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # proper_motion_component_vr
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # right_ascension
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # declination
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # galactocentric_distance
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # galactocentric_coordinate_b
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # galactocentric_coordinate_l
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # go_photometry
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # gr_photometry
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # rz_photometry
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # v_photometry
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # velocity_u
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # velocity_v
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # velocity_w
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # spectral_type
    strategies.integers())

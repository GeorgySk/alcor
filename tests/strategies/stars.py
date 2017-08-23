import sys
from decimal import Decimal

from hypothesis import strategies

from alcor.models import Star
from alcor.models.star import GalacticDiskEnum
from alcor.services.common import STARS_SPECTRAL_TYPES

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
    # galactic_distance
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # galactic_latitude
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # galactic_longitude
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # ugriz_g_apparent
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # ugriz_ug
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # ugriz_gr
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # ugriz_ri
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # ugriz_iz
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # v_photometry
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # u_velocity
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # v_velocity
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # w_velocity
    strategies.decimals(allow_nan=False,
                        allow_infinity=False),
    # spectral_type
    strategies.one_of(*map(strategies.just, STARS_SPECTRAL_TYPES)),
    # disk_belonging
    strategies.one_of(*map(strategies.just, GalacticDiskEnum)))

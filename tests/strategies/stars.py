import sys
from decimal import Decimal

from hypothesis import strategies
from hypothesis.searchstrategy.strategies import MappedSearchStrategy

from alcor.models import Star
from alcor.models.star import GalacticDiskType
from alcor.services.common import STARS_SPECTRAL_TYPES


def stars_factory(nullable: bool) -> MappedSearchStrategy:
    decimals = strategies.decimals(allow_nan=False,
                                   allow_infinity=False)
    positive_decimals = strategies.decimals(
        min_value=Decimal(sys.float_info.epsilon),
        allow_nan=False,
        allow_infinity=False)

    if nullable:
        decimals |= strategies.none()
        positive_decimals |= strategies.none()

    return strategies.builds(
        Star,
        # group_id
        strategies.text(min_size=6),
        # luminosity
        decimals,
        # proper_motion
        positive_decimals,
        # proper_motion_component_b
        decimals,
        # proper_motion_component_l
        decimals,
        # proper_motion_component_vr
        decimals,
        # right_ascension
        decimals,
        # declination
        decimals,
        # galactic_distance
        decimals,
        # galactic_latitude
        decimals,
        # galactic_longitude
        decimals,
        # ugriz_g_apparent
        decimals,
        # ugriz_ug
        decimals,
        # ugriz_gr
        decimals,
        # ugriz_ri
        decimals,
        # ugriz_iz
        decimals,
        # v_photometry
        decimals,
        # u_velocity
        decimals,
        # v_velocity
        decimals,
        # w_velocity
        decimals,
        # spectral_type
        strategies.one_of(map(strategies.just, STARS_SPECTRAL_TYPES)),
        # galactic_disk_type
        strategies.one_of(map(strategies.just, GalacticDiskType)))


defined_stars = stars_factory(nullable=False)
undefined_stars = stars_factory(nullable=True)

defined_stars_lists = strategies.lists(defined_stars)
undefined_stars_lists = strategies.lists(undefined_stars)

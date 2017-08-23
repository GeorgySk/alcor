from collections import Counter
from decimal import Decimal
from fractions import Fraction
from math import log10

from alcor.models import Star


def check(star: Star,
          *,
          # we are using `dict`s mutability,
          # so passed as argument
          # original object will have all changes
          # which took place inside of a function
          eliminations_counter: Counter,
          filtration_method: str,
          min_parallax: float = 0.025,
          min_declination: float = 0.,
          max_velocity: float = 500.,
          min_proper_motion: float = 0.04) -> bool:
    # TODO: implement pc/kpc units
    galactic_distance = star.galactic_distance * Decimal(1e3)
    parallax = Fraction(1, Fraction(galactic_distance))
    # TODO: find out the meaning of the following constants
    hrm = star.ugriz_g_apparent + Decimal(5. * log10(star.proper_motion) + 5.)
    gz = float(star.ugriz_gr) + float(star.ugriz_rz)

    if parallax < min_parallax:
        eliminations_counter['parallax'] += 1
        return True

    northern_hemisphere_star = star.declination < min_declination
    if northern_hemisphere_star:
        eliminations_counter['declination'] += 1
        return True

    hypervelocity_star = (star.u_velocity ** 2
                          + star.v_velocity ** 2
                          + star.w_velocity ** 2
                          > max_velocity ** 2)
    if hypervelocity_star:
        eliminations_counter['velocity'] += 1
        return True

    if filtration_method == 'restricted':
        if star.proper_motion < min_proper_motion:
            eliminations_counter['proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif gz < -0.33 and hrm < 14.:
            eliminations_counter['reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif hrm < 3.559 * gz + 15.17:
            eliminations_counter['reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constant
        elif star.v_photometry >= 19.:
            eliminations_counter['apparent_magnitude'] += 1
            return True
    return False

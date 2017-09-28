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
    distance_in_pc = star.distance * Decimal(1e3)
    parallax = Fraction(1, Fraction(distance_in_pc))

    if parallax < min_parallax:
        eliminations_counter['by_parallax'] += 1
        return True

    northern_hemisphere_star = star.declination < min_declination
    if northern_hemisphere_star:
        eliminations_counter['by_declination'] += 1
        return True

    hypervelocity_star = (star.u_velocity ** 2
                          + star.v_velocity ** 2
                          + star.w_velocity ** 2
                          > max_velocity ** 2)
    if hypervelocity_star:
        eliminations_counter['by_velocity'] += 1
        return True

    if filtration_method == 'restricted':
        # TODO: add properties or function for converting?
        # Transformation from UBVRI to ugriz. More info at:
        # Jordi, Grebel & Ammon, 2006, A&A, 460; equations 1-8 and Table 3
        g_ugriz_abs_magnitude = (float(star.v_abs_magnitude) - 0.124
                                 + 0.63 * float(star.b_abs_magnitude
                                                - star.v_abs_magnitude))
        z_ugriz_abs_magnitude = (g_ugriz_abs_magnitude
                                 - 1.646 * float(star.v_abs_magnitude
                                                 - star.r_abs_magnitude)
                                 - 1.584 * float(star.r_abs_magnitude
                                                 - star.i_abs_magnitude)
                                 + 0.525)
        g_apparent_magnitude = apparent_magnitude(
                g_ugriz_abs_magnitude,
                distance_kpc=float(star.distance))
        z_apparent_magnitude = apparent_magnitude(
                z_ugriz_abs_magnitude,
                distance_kpc=float(star.distance))
        # TODO: find out the meaning and check if the last 5 is correct
        hrm = g_apparent_magnitude + float(
                5. * log10(star.proper_motion) + 5.)
        v_apparent_magnitude = apparent_magnitude(
                float(star.v_abs_magnitude),
                distance_kpc=float(star.distance))

        if star.proper_motion < min_proper_motion:
            eliminations_counter['by_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif (g_apparent_magnitude - z_apparent_magnitude < -0.33
              and hrm < 14.):
            eliminations_counter['by_reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif (hrm < 3.559 * (g_apparent_magnitude - z_apparent_magnitude)
              + 15.17):
            eliminations_counter['by_reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constant
        elif v_apparent_magnitude >= 19.:
            eliminations_counter['by_apparent_magnitude'] += 1
            return True
    return False


def apparent_magnitude(abs_magnitude: float,
                       distance_kpc: float
                       ) -> float:
    # More info at (2nd formula, + 3.0 because the distance is in kpc):
    # https://en.wikipedia.org/wiki/Absolute_magnitude#Apparent_magnitude
    return abs_magnitude - 5. + 5. * (log10(distance_kpc) + 3.)

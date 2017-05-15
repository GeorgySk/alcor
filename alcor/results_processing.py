from collections import Counter
from fractions import Fraction
from functools import partial
from itertools import filterfalse
from math import (log10)
from typing import Dict

import astropy.units as u

from alcor.models.star import Star
from alcor.utils import parse_stars


def main() -> None:

    # Binning parameters of Luminosity Function
    min_bolometric_magnitude = 6.0
    max_bolometric_magnitude = 21.0
    bin_size = 0.5
    bolometric_magnitude_amplitude = (max_bolometric_magnitude
                                      - min_bolometric_magnitude)
    bins_count = bolometric_magnitude_amplitude / bin_size

    with open('output.res', 'r') as output_file:
        full_stars_sample = list(parse_stars(output_file, 'test'))

    elimination_counters = Counter(int)
    apply_elimination_criteria = partial(check_if_eliminated,
                                         elimination_counters
                                         =elimination_counters)
    restricted_stars_sample = filterfalse(apply_elimination_criteria,
                                          full_stars_sample)

    for star in restricted_stars_sample:
        star.set_radial_velocity_to_zero()
        distribute_into_bins(star, bins)


def check_if_eliminated(star: Star,
                        elimination_counters: Dict[str, int]) -> bool:
    min_parallax = 0.025 * u.arcsec
    min_declination = 0.0
    max_velocity = 500.0
    min_proper_motion = 0.04

    galactocentric_distance = star.galactocentric_distance * 10e3  # From kpc
    parallax = Fraction(1.0, Fraction(galactocentric_distance))
    hrm = star.go_photometry + 5.0 * log10(star.proper_motion) + 5.0
    gz = star.gr_photometry + star.rz_photometry

    if parallax < min_parallax:
        elimination_counters['parallax'] += 1
        return True
    elif star.declination < min_declination:
        elimination_counters['declination'] += 1
        return True
    elif (star.velocity_u**2 + star.velocity_v**2 + star.velocity_z**2
            > max_velocity**2):
        elimination_counters['velocity'] += 1
        return True
    elif star.proper_motion < min_proper_motion:
        elimination_counters['proper_motion'] += 1
        return True
    elif gz < -0.33 and hrm < 14.0:
        elimination_counters['reduced_proper_motion'] += 1
        return True
    elif hrm < 3.559 * gz + 15.17:
        elimination_counters['reduced_proper_motion'] += 1
        return True
    elif star.v_photometry >= 19.0:
        elimination_counters['apparent_magnitude'] += 1
        return True
    return False


if __name__ == '__main__':
    main()

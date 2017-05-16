from collections import Counter
from fractions import Fraction
from typing import Dict

from math import log10

from alcor.models.star import Star

MIN_PARALLAX = 0.025
MIN_DECLINATION = 0.0
MAX_VELOCITY = 500.0
MIN_PROPER_MOTION = 0.04


def check_if_eliminated(star: Star,
                        elimination_counters: Dict[str, int]) -> bool:
    galactocentric_distance = star.galactocentric_distance * 10e3
    parallax = Fraction(1.0, Fraction(galactocentric_distance))
    hrm = star.go_photometry + 5.0 * log10(star.proper_motion) + 5.0
    gz = star.gr_photometry + star.rz_photometry

    if parallax < MIN_PARALLAX:
        elimination_counters['parallax'] += 1
        return True
    elif star.declination < MIN_DECLINATION:
        elimination_counters['declination'] += 1
        return True
    elif (star.velocity_u ** 2 + star.velocity_v ** 2 + star.velocity_z ** 2
          > MAX_VELOCITY ** 2):
        elimination_counters['velocity'] += 1
        return True
    elif star.proper_motion < MIN_PROPER_MOTION:
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


def write_elimination_stats(full_sample_stars_count: int,
                            restricted_sample_stars_count: int,
                            elimination_counters: Counter[int]) -> None:
    with open('elimination_stats.res', 'w') as output_file:
        output_file.write('Initial number of WDs:',
                          full_sample_stars_count,
                          'Eliminated by parallax:',
                          elimination_counters[parallax],
                          'Eliminated by declination:',
                          elimination_counters[declination],
                          'Initial number of stars in northern hemisphere:',
                          full_sample_stars_count
                          - elimination_counters[parallax]
                          - elimination_counters[declination],
                          'Eliminated by proper motion:',
                          elimination_counters[proper_motion],
                          'Eliminated by reduced proper motion:',
                          elimination_counters[reduced_proper_motion],
                          'Eliminated by apparent magnitude:',
                          elimination_counters[apparent_magnitude],
                          'Number of stars in restricted sample:',
                          restricted_sample_stars_count)
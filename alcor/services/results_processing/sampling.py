import csv
from collections import Counter
from fractions import Fraction

from math import log10

from alcor.models.star import Star

MIN_PARALLAX = 0.025
MIN_DECLINATION = 0.0
MAX_VELOCITY = 500.0
MIN_PROPER_MOTION = 0.04


def check_elimination(star: Star,
                      elimination_counters: Counter,
                      method: str) -> bool:
    # TODO: implement pc/kpc units
    galactocentric_distance = star.galactocentric_distance * 10e3
    parallax = Fraction(1.0, Fraction(galactocentric_distance))
    # TODO: find out the meaning of the following constants
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
    elif method == 'restricted':
        if star.proper_motion < MIN_PROPER_MOTION:
            elimination_counters['proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif gz < -0.33 and hrm < 14.0:
            elimination_counters['reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif hrm < 3.559 * gz + 15.17:
            elimination_counters['reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constant
        elif star.v_photometry >= 19.0:
            elimination_counters['apparent_magnitude'] += 1
            return True
    return False


def write_elimination_stats(full_sample_stars_count: int,
                            restricted_sample_stars_count: int,
                            elimination_counters: Counter) -> None:
    with open('elimination_stats.csv', 'w') as csv_file:
        file_writer = csv.writer(csv_file, delimiter='  ')
        file_writer.writerow('Initial number of WDs:',
                             'Eliminated by parallax:',
                             'Eliminated by declination:',
                             'Initial number of stars in northern hemisphere:',
                             'Eliminated by proper motion:',
                             'Eliminated by reduced proper motion:',
                             'Eliminated by apparent magnitude:',
                             'Number of stars in restricted sample:')
        file_writer.writerow(full_sample_stars_count,
                             elimination_counters['parallax'],
                             elimination_counters['declination'],
                             full_sample_stars_count
                             - elimination_counters['parallax']
                             - elimination_counters['declination'],
                             elimination_counters['proper_motion'],
                             elimination_counters['reduced_proper_motion'],
                             elimination_counters['apparent_magnitude'],
                             restricted_sample_stars_count)

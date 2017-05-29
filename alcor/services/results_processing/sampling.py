import csv
import logging
from collections import Counter
from fractions import Fraction

from math import log10

from cassandra.deserializers import Decimal

from alcor.models.star import Star

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)

MIN_PARALLAX = 0.025
MIN_DECLINATION = 0.
MAX_VELOCITY = 500.
MIN_PROPER_MOTION = 0.04


def check_elimination(star: Star,
                      eliminations_counter: Counter,
                      method: str) -> bool:
    # TODO: implement pc/kpc units
    galactocentric_distance = star.galactocentric_distance * Decimal(1e3)
    parallax = Fraction(1, Fraction(galactocentric_distance))
    # TODO: find out the meaning of the following constants
    hrm = star.go_photometry + Decimal(5. * log10(star.proper_motion) + 5.)
    gz = star.gr_photometry + star.rz_photometry

    if parallax < MIN_PARALLAX:
        eliminations_counter['parallax'] += 1
        return True
    elif star.declination < MIN_DECLINATION:
        eliminations_counter['declination'] += 1
        return True
    elif (star.velocity_u ** 2 + star.velocity_v ** 2 + star.velocity_w ** 2
          > MAX_VELOCITY ** 2):
        eliminations_counter['velocity'] += 1
        return True
    elif method == 'restricted':
        if star.proper_motion < MIN_PROPER_MOTION:
            eliminations_counter['proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif gz < -0.33 and hrm < 14.:
            eliminations_counter['reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constants
        elif hrm < Decimal(3.559) * gz + Decimal(15.17):
            eliminations_counter['reduced_proper_motion'] += 1
            return True
        # TODO: find out the meaning of the following constant
        elif star.v_photometry >= 19.:
            eliminations_counter['apparent_magnitude'] += 1
            return True
    return False


def write_elimination_stats(raw_sample_stars_count: int,
                            filtered_stars_count: int,
                            eliminations_counter: Counter) -> None:
    with open('elimination_stats.csv', 'w') as csv_file:
        file_writer = csv.writer(csv_file, delimiter=' ')
        file_writer.writerow(['RawNº',
                              'Par-x',
                              'Dec',
                              'NH Nº',
                              'PropM',
                              'RedPrM',
                              'AppMag',
                              'RestrNº'])
        file_writer.writerow([raw_sample_stars_count,
                              eliminations_counter['parallax'],
                              eliminations_counter['declination'],
                              raw_sample_stars_count
                              - eliminations_counter['parallax']
                              - eliminations_counter['declination'],
                              eliminations_counter['proper_motion'],
                              eliminations_counter['reduced_proper_motion'],
                              eliminations_counter['apparent_magnitude'],
                              filtered_stars_count])

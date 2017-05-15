from collections import defaultdict
from fractions import Fraction
from math import (log10,
                  sqrt,
                  cos,
                  sin)
import quantities as pq
from typing import Dict

from alcor.models.star import Star
from alcor.utils import parse_stars


def main() -> None:
    stars = []

    # Binning parameters of Luminosity Function
    min_bolometric_magnitude = 6.0
    max_bolometric_magnitude = 21.0
    bin_size = 0.5
    bolometric_magnitude_amplitude = (max_bolometric_magnitude
                                      - min_bolometric_magnitude)
    bins_count = bolometric_magnitude_amplitude / bin_size

    with open('output.res', 'r') as output_file:
        stars = list(parse_stars(output_file, 'test'))

    elimination_counters = defaultdict(int)

    for star in stars:
        if not is_eliminated(star,
                             elimination_counters):
            star.set_radial_velocity_to_zero()
            distribute_into_bins(star, bins)


def is_eliminated(star: Star,
                  elimination_counters: Dict[str, int]) -> bool:
    min_parallax = 0.025
    min_declination = 0.0
    max_velocity = 500.0
    min_proper_motion = 0.04

    galactocentric_distance = star.galactocentric_distance * 10e3  # From kpc
    parallax = Fraction(1.0, galactocentric_distance)
    hrm = star.go_photometry + 5.0*log10(star.proper_motion) + 5.0
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
    elif (gz < -0.33) and (hrm < 14.0):
        elimination_counters['reduced_proper_motion'] += 1
        return True
    elif hrm < (3.559 * gz + 15.17):
        elimination_counters['reduced_proper_motion'] += 1
        return True
    elif star.v_photometry >= 19.0:
        elimination_counters['apparent_magnitude'] += 1
        return True
    return False


if __name__ == '__main__':
    main()

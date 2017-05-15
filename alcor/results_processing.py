from collections import defaultdict
from decimal import Decimal
from math import (log10,
                  sqrt)
from typing import Dict

from alcor.models.star import Star


def main() -> None:
    stars = []

    # Binning parameters of Luminosity Function
    min_bolometric_magnitude = 6.0
    max_bolometric_magnitude = 21.0
    bin_size = 0.5
    bolometric_magnitude_amplitude = (max_bolometric_magnitude
                                      - min_bolometric_magnitude)
    bins_count = bolometric_magnitude_amplitude / bin_size

    with open('output.res', 'r') as f:
        line_counter = 0
        for line in f:
            line_counter += 1
            parts = line.split()
            params = map(Decimal, parts)
            star = Star(*params)
            stars.append(star)

    elimination_counters = defaultdict(int)

    for star in stars:
        star_is_eliminated = False
        apply_elimination_criteria(star,
                                   star_is_eliminated,
                                   elimination_counters)
        if not star_is_eliminated:
            # Save stars after applying elimination criteria
            set_radial_velocity_to_zero(stars)
            distribute_into_bins(star, bins)


def apply_elimination_criteria(star: Star,
                               star_is_eliminated: bool,
                               elimination_counters: Dict[str, int]):
    min_parallax = 0.025
    min_declination = 0.0
    max_velocity = 500.0
    min_proper_motion = 0.04

    galactocentric_distance = star.galactocentric_distance * 10e3  # From kpc
    parallax = 1.0 / galactocentric_distance
    hrm = star.go_photometry + 5.0*log10(star.proper_motion) + 5.0
    gz = star.gr_photometry + star.rz_photometry

    if parallax < min_parallax:
        elimination_counters["eliminated_by_parallax"] += 1
        star_is_eliminated = True
    elif star.declination < min_declination:
        elimination_counters["eliminated_by_declination"] += 1
        star_is_eliminated = True
    elif (sqrt(star.velocity_u**2 + star.velocity_v**2 + star.velocity_z**2)
            > max_velocity):
        elimination_counters["eliminated_by_velocity"] += 1
        star_is_eliminated = True
    elif star.proper_motion < min_proper_motion:
        elimination_counters["eliminated_by_proper_motion"] += 1
        star_is_eliminated = True
    elif (gz < -0.33) and (hrm < 14.0):
        elimination_counters["eliminated_by_reduced_proper_motion"] += 1
        star_is_eliminated = True
    elif hrm < (3.559 * gz + 15.17):
        elimination_counters["eliminated_by_reduced_proper_motion"] += 1
        star_is_eliminated = True
    elif star.v_photometry >= 19.0:
        elimination_counters["eliminated_by_apparent_magnitude"] += 1
        star_is_eliminated = True
    return elimination_counters, star_is_eliminated


if __name__ == '__main__':
    main()

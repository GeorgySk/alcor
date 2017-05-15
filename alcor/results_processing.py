from collections import namedtuple
from decimal import Decimal
from math import log10, sqrt

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

    elimination_counters = {'eliminated_by_parallax': 0,
                            'eliminated_by_declination': 0,
                            'eliminated_by_velocity': 0,
                            'eliminated_by_proper_motion': 0,
                            'eliminated_by_reduced_proper_motion': 0,
                            'eliminated_by_apparent_magnitude': 0}

    for i in range(0, size(stars)):
        star_is_eliminated = False
        apply_elimination_criteria(stars(i),
                                   star_is_eliminated,
                                   elimination_counters)
        if not star_is_eleminated:
            # Save stars after applying elimination criteria
            set_radial_velocity_to_zero(stars(i))
            distribute_into_bins(star, bins)


def apply_elimination_criteria(star: Star,
                               star_is_eliminated: bool,
                               elimination_counters: dict):
    min_parallax = 0.025
    min_declination = 0.0
    max_velocity = 500.0
    min_proper_motion = 0.04

    galactocentric_distance = star.galactocentric_distance * 10e3  # From kpc
    parallax = 1.0 / galactocentric_distance
    tangential_velocity = 4.74 * star.proper_motion * galactocentric_distance
    hrm = star.go_photometry + 5.0*log10(star.proper_motion) + 5.0
    gz = star.gr_photometry + star.rz_photometry

    if  parallax < min_parallax:
        elimination_counters(eliminated_by_parallax) += 1
        star_is_eliminated = True
    else if star.declination < min_declination:
        elimination_counters(eliminated_by_declination) += 1
        star_is_eliminated = True
    else if sqrt(star.velocity_u**2 + star.velocity_v**2 + star.velocity_z**2) > max_velocity:
        elimination_counters(eliminated_by_velocity) += 1
        star_is_eliminated = True
    else if star.proper_motion < min_proper_motion:
        elemination_counters(eliminated_by_proper_motion) += 1
        star_is_eliminated = True
    else if (gz < -0.33) and (hrm < 14.0):
        elemination_counters(eliminated_by_reduced_proper_motion) += 1
        star_is_eliminated = True
    else if hrm < (3.559 * gz + 15.17):
        elemination_counters(eliminated_by_reduced_proper_motion) += 1
        star_is_eliminated = True
    else if star.v_photometry >= 19.0:
        elemination_counters(eliminated_by_apparent_magnitude) += 1
        star_is_eliminated = True


if __name__ == '__main__':
    main()

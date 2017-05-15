from math import atan
from collections import namedtuple
from decimal import Decimal


def main() -> None:
    stars = []
    Star = namedtuple('Star', 'luminosity '
                              'proper_motion '
                              'declination ' 
                              'galactocentric_distance '
                              'gr_photometry '
                              'rz_photometry '
                              'v_photometry '
                              'velocity_u '
                              'velocity_v '
                              'velocity_w ')

    # Binning parameters of Luminosity Function
    min_bolometric_magnitude = 6.0
    max_bolometric_magnitude = 21.0
    bin_size = 0.5
    bolometric_magnitude_amplitude = (max_bolometric_magnitude
                                      - min_bolometric_magnitude)
    bins_count = bolometric_magnitude_amplitude / bin_size

    with open('output.dat', 'r') as f:
        line_counter = 0
        for line in f:
            line_counter += 1
            parts = line.split()
            params = map(Decimal, parts)
            star = Star(*params)
            stars.append(star)

    elemination_counters = {'eliminated_by_parallax': 0,
                            'eliminated_by_declination': 0,
                            'eliminated_by_velocity': 0,
                            'eliminated_by_proper_motion': 0,
                            'eliminated_by_reduced_proper_motion': 0,
                            'eliminated_by_apparent_magnitude': 0}

    for i in range(0, size(stars)):
        star_is_eleminated = False
        apply_elemination_criteria(stars(i), star_is_eleminated, elemination_counters)
        if not star_is_eleminated:
            # Save stars after applying elemination criteria
            set_radial_velocity_to_zero(stars(i))
            distribute_into_bins(star, bins)

            



if __name__ == '__main__':
    main()

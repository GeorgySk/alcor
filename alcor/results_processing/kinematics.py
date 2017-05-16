from alcor.models.star import Star
from alcor.results_processing.luminosity_function import \
    (MIN_BOLOMETRIC_MAGNITUDE,
     BIN_SIZE)


def write_bins_kinematic_info(bins: list[list[Star]]) -> None:
    with open('magnitude_bins.res', 'w') as output_file:
        for bin_index, bin in enumerate(bins):
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (bin_index - 0.5))
            average_bin_velocity_u = sum(bin(:).velocity_u) / len(bin)
            average_bin_velocity_v = sum(bin(:).velocity_v) / len(bin)
            average_bin_velocity_w = sum(bin(:).velocity_w) / len(bin)
            bin_standard_deviation_of_velocity_u = std.(bin(:).velocity_u)
            bin_standard_deviation_of_velocity_v = std.(bin(:).velocity_v)
            bin_standard_deviation_of_velocity_w = std.(bin(:).velocity_w)
            output_file.write(average_bin_magnitude,
                              average_bin_velocity_u,
                              average_bin_velocity_v,
                              average_bin_velocity_w,
                              bin_standard_deviation_of_velocity_u,
                              bin_standard_deviation_of_velocity_v,
                              bin_standard_deviation_of_velocity_w)
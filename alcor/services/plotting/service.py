import csv
from operator import (add,
                      sub)

from bokeh.plotting import (figure,
                            output_file,
                            show)


def make_plots() -> None:
    with open('luminosity_function.csv', 'r') as data_file:
        reader = csv.reader(data_file)

        first_row = next(reader)
        splitted_first_row = first_row[0].split()
        normalization_factor = float(splitted_first_row[1])

        header_row = next(reader)

        average_bin_magnitude = []
        star_count_logarithm = []
        upper_errorbar = []
        lower_errorbar = []

        for row in reader:
            splitted_row = row[0].split()
            average_bin_magnitude.append(float(splitted_row[0]))
            star_count_logarithm.append(float(splitted_row[1]))
            upper_errorbar.append(float(splitted_row[2]))
            lower_errorbar.append(float(splitted_row[3]))

            output_file("luminosity_function.html")

        p = figure(title="luminosity function",
                   x_axis_label='Mbol',  # Todo: rename axis
                   y_axis_label='logN')

        # TODO: find out if we can omit legend here
        p.line(average_bin_magnitude,
               star_count_logarithm,
               legend="Temp.",
               line_width=2)
        p.square(average_bin_magnitude,
                 star_count_logarithm)
        p.multi_line(average_bin_magnitude,
                     map(add,
                         star_count_logarithm,
                         upper_errorbar))
        p.multi_line(average_bin_magnitude,
                     map(sub,
                         star_count_logarithm,
                         lower_errorbar))

        show(p)

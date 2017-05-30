import csv

from bokeh.plotting import (figure,
                            output_file,
                            show)


def plot_luminosity_function() -> None:
    output_file("luminosity_function.html")

    plot = figure(title="luminosity function",
                  x_axis_label='Mbol',  # Todo: rename axis
                  y_axis_label='logN')

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

    # TODO: find out if we can omit legend here
    plot.line(average_bin_magnitude,
              star_count_logarithm,
              legend="Temp.",
              line_width=2)
    plot.square(average_bin_magnitude,
                star_count_logarithm)

    multiline_x = []
    multiline_y = []

    for (magnitude, stars_count,
         error_up, error_down) in zip(average_bin_magnitude,
                                      star_count_logarithm,
                                      upper_errorbar,
                                      lower_errorbar):
        multiline_x.append((magnitude,
                            magnitude))
        multiline_y.append((stars_count - error_down,
                            stars_count + error_up))

    plot.multi_line(multiline_x,
                    multiline_y)

    show(plot)

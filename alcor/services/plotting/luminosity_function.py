import csv

from bokeh.plotting import (figure,
                            output_file,
                            show)


def plot() -> None:
    output_file("luminosity_function.html")

    main_plot = figure(title="luminosity function",
                       x_axis_label='Mbol',  # Todo: rename axis
                       y_axis_label='logN')

    with open('luminosity_function.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=' ')

        first_row = next(reader)
        normalization_factor = float(first_row[1])

        header_row = next(reader)

        average_bin_magnitude = []
        star_count_logarithm = []
        upper_errorbar = []
        lower_errorbar = []

        for row in reader:
            average_bin_magnitude.append(float(row[0]))
            star_count_logarithm.append(float(row[1]))
            upper_errorbar.append(float(row[2]))
            lower_errorbar.append(float(row[3]))

    # TODO: find out if we can omit legend here
    main_plot.line(average_bin_magnitude,
                   star_count_logarithm,
                   legend="Temp.",
                   line_width=2)
    main_plot.square(average_bin_magnitude,
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

    main_plot.multi_line(multiline_x,
                         multiline_y)

    show(main_plot)

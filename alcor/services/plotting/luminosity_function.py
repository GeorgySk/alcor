import csv
from math import pi
from typing import List

from bokeh.plotting.figure import Figure

from alcor.utils import get_columns
from .latex_label import LatexLabel

from bokeh.plotting import (figure,
                            output_file,
                            show)

CSV_DELIMITER = ' '


def plot() -> None:
    output_file("luminosity_function.html")

    main_plot = figure(title="Luminosity function")

    # TODO: get coordinates for labels automatically
    add_latex_labels(fig=main_plot,
                     x_label_text='M_{bol}',
                     x_label_xpos=300,
                     x_label_ypos=0,
                     y_label_text='\log N (pc^{-3}M_{bol}^{-1})',
                     y_label_xpos=-60,
                     y_label_ypos=300)

    with open('luminosity_function.csv', 'r') as file:
        reader = csv.reader(file,
                            delimiter=CSV_DELIMITER)
        first_row = next(reader)
        # TODO: figure out what do do with normalization_factor
        normalization_factor = float(first_row[1])
        header_row = next(reader)

        (average_bin_magnitude,
         star_count_logarithm,
         upper_errorbar,
         lower_errorbar) = get_columns(reader)

    main_plot.line(average_bin_magnitude,
                   star_count_logarithm,
                   line_width=2)
    main_plot.square(average_bin_magnitude,
                     star_count_logarithm)

    add_errorbars(fig=main_plot,
                  x=average_bin_magnitude,
                  y=star_count_logarithm,
                  y_err_up=upper_errorbar,
                  y_err_down=lower_errorbar)

    show(main_plot)


def add_latex_labels(fig: Figure,
                     x_label_text: str,
                     x_label_xpos: int,
                     x_label_ypos: int,
                     y_label_text: str,
                     y_label_xpos: int,
                     y_label_ypos: int) -> None:
    x_label_latex = LatexLabel(text=x_label_text,
                               x=x_label_xpos,
                               y=x_label_ypos,
                               x_units='screen',
                               y_units='screen',
                               render_mode='css',
                               text_font_size='8pt',
                               background_fill_color='#ffffff')
    fig.add_layout(x_label_latex)
    y_label_latex = LatexLabel(text=y_label_text,
                               x=y_label_xpos,
                               y=y_label_ypos,
                               angle=pi / 2.,
                               x_units='screen',
                               y_units='screen',
                               render_mode='css',
                               text_font_size='8pt',
                               background_fill_color='#ffffff')
    fig.add_layout(y_label_latex)


def add_errorbars(fig: Figure,
                  x: List[float],
                  y: List[float],
                  y_err_up: List[float],
                  y_err_down: List[float]) -> None:
    multiline_x = []
    multiline_y = []
    for (x_coordinate, y_coordinate,
         error_up, error_down) in zip(x, y, y_err_up, y_err_down):
        multiline_x.append((x_coordinate,
                            x_coordinate))
        multiline_y.append((y_coordinate - error_down,
                            y_coordinate + error_up))

    fig.multi_line(multiline_x,
                   multiline_y)

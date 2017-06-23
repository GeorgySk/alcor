import csv
from math import pi
from typing import List

from bokeh.plotting.figure import Figure
from cassandra.cluster import Session

from alcor.models.luminosity_function import Point
from alcor.services.data_access.reading import fetch
from alcor.utils import get_columns
from .latex_label import LatexLabel

from bokeh.plotting import (figure,
                            output_file,
                            show)

CSV_DELIMITER = ' '


def plot(*,
         session: Session) -> None:

    graph_points = fetch_last_graph_points(session=session)

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

    main_plot.line(graph_points.avg_bin_magnitude,
                   graph_points.star_count_logarithm,
                   line_width=2)
    main_plot.square(graph_points.avg_bin_magnitude,
                     graph_points.star_count_logarithm)

    add_errorbars(fig=main_plot,
                  x=graph_points.average_bin_magnitude,
                  y=graph_points.star_count_logarithm,
                  y_err_up=graph_points.upper_error_bar,
                  y_err_down=graph_points.lower_error_bar)

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


def fetch_last_graph_points(*,
                            session: Session) -> List[Point]:
    query = (Point.objects.all().order_by('-updated_timestamp'))
    query = query[0]
    records = fetch(query=query,
                    session=session)
    return [Point(**record)
            for record in records]

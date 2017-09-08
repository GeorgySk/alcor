import logging
from typing import (Tuple,
                    List)

from matplotlib.axes import Axes
from sqlalchemy.orm.session import Session
import matplotlib
import numpy as np

from alcor.models.star import Star
from alcor.services.data_access import (fetch_all,
                                        fetch_all_da_stars)
# More info at
#  http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot(session: Session,
         ug_gr_filename: str = 'ug_gr_ugriz.png',
         gr_ri_filename: str = 'gr_ri_ugriz.png',
         ri_iz_filename: str = 'ri_iz_ugriz.png',
         gr_uggr_limit: Tuple[float, float] = (-0.8, 0.5),
         ug_uggr_limit: Tuple[float, float] = (1., -0.8),
         ri_grri_limit: Tuple[float, float] = (-0.65, 0.2),
         gr_grri_limit: Tuple[float, float] = (0.5, -0.8),
         ri_riiz_limit: Tuple[float, float] = (0.5, -0.8),
         iz_riiz_limit: Tuple[float, float] = (-1.7, 1.1),
         ug_label: str = '$u-g$',
         gr_label: str = '$g-r$',
         ri_label: str = '$r-i$',
         iz_label: str = '$i-z$',
         da_only: bool = True) -> None:
    # TODO: add more fetching options
    if da_only:
        stars = fetch_all_da_stars(session)
    else:
        stars = fetch_all(Star,
                          session=session)

    ugriz_ug = [float(star.ugriz_ug)
                for star in stars]
    ugriz_gr = [float(star.ugriz_gr)
                for star in stars]
    ugriz_ri = [float(star.ugriz_ri)
                for star in stars]
    ugriz_iz = [float(star.ugriz_iz)
                for star in stars]

    draw_plot(xlabel=gr_label,
              ylabel=ug_label,
              x=ugriz_gr,
              y=ugriz_ug,
              xlim=gr_uggr_limit,
              ylim=ug_uggr_limit,
              filename=ug_gr_filename,
              constraints_contours='ug-gr')
    draw_plot(xlabel=ri_label,
              ylabel=gr_label,
              x=ugriz_ri,
              y=ugriz_gr,
              xlim=ri_grri_limit,
              ylim=gr_grri_limit,
              filename=gr_ri_filename,
              constraints_contours='gr-ri')
    draw_plot(xlabel=iz_label,
              ylabel=ri_label,
              x=ugriz_iz,
              y=ugriz_ri,
              xlim=iz_riiz_limit,
              ylim=ri_riiz_limit,
              filename=ri_iz_filename,
              constraints_contours='ri-iz')


def draw_plot(xlabel: str,
              ylabel: str,
              x: List[float],
              y: List[float],
              xlim: Tuple[float, float],
              ylim: Tuple[float, float],
              filename: str,
              constraints_contours: str,
              figure_size: Tuple[float, float] = (8, 8),
              color: str = 'b',
              point_size: float = 0.5,
              ratio: float = 10 / 13
              ) -> None:
    figure, subplot = plt.subplots(figsize=figure_size)

    subplot.set(xlabel=xlabel,
                ylabel=ylabel,
                xlim=xlim,
                ylim=ylim)

    subplot.scatter(x=x,
                    y=y,
                    color=color,
                    s=point_size)

    if constraints_contours == 'ug-gr':
        add_ug_gr_constraints_contours(subplot=subplot)
    elif constraints_contours == 'gr-ri':
        add_gr_ri_constraints_contours(subplot=subplot)
    elif constraints_contours == 'ri-iz':
        add_ri_iz_constraints_contours(subplot=subplot)

    subplot.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    plt.savefig(filename)


# More info on polynomials at 'DA white dwarfs in Sloan Digital Sky Survey Data
# Release 7 and a search for infrared excess emission. J.Girven et al.'
def add_ug_gr_constraints_contours(subplot: Axes) -> None:
    upper_polynomial = np.poly1d([-24.384, -19.,
                                  3.497, 1.193,
                                  0.083, 0.61])
    lower_polynomial = np.poly1d([-20.653, 10.816,
                                  15.718, -1.294,
                                  -0.084, 0.3])
    # Approx. intersection locations can be found from Fig.2 of Girven's
    # article. Exact values from:
    # (upper_polynomial - lower_polynomial).roots
    gr = np.linspace(-0.62128319, 0.30067351, 100)

    ug_upper = upper_polynomial(gr)
    ug_lower = lower_polynomial(gr)

    subplot.plot(gr, ug_upper,
                 color='r')
    subplot.plot(gr, ug_lower,
                 color='r')


def add_gr_ri_constraints_contours(subplot: Axes) -> None:
    upper_polynomial = np.poly1d([-0.6993, 0.947, 0.192])
    lower_polynomial = np.poly1d([-1.32, 2.173, 2.452, -0.07])

    # Approx. intersection locations can be found from Fig.2 of Girven's
    # article. Exact values from:
    # (upper_polynomial - lower_polynomial).roots
    ri = np.linspace(-0.55044287, 0.13938346, 100)

    gr_upper = upper_polynomial(ri)
    gr_lower = lower_polynomial(ri)

    subplot.plot(ri, gr_upper,
                 color='r')
    subplot.plot(ri, gr_lower,
                 color='r')


def add_ri_iz_constraints_contours(subplot: Axes) -> None:
    top_line = np.poly1d([-0.56])
    left_line = np.poly1d([0.176, 0.127])
    right_line = np.poly1d([-0.754, 0.11])

    right_nod = (top_line - right_line).roots[0]
    middle_nod = (left_line - right_line).roots[0]
    left_nod = (top_line - left_line).roots[0]

    iz_upper_line = [left_nod, right_nod]
    iz_left_line = [left_nod, middle_nod]
    iz_right_line = [middle_nod, right_nod]

    ri_upper = top_line(iz_upper_line)
    ri_left = left_line(iz_left_line)
    ri_right = right_line(iz_right_line)

    subplot.plot(iz_upper_line, ri_upper,
                 color='r')
    subplot.plot(iz_left_line, ri_left,
                 color='r')
    subplot.plot(iz_right_line, ri_right,
                 color='r')

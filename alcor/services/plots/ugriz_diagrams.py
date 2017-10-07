import logging
from typing import (Tuple,
                    List)

import matplotlib
from matplotlib.axes import Axes

# More info at
#  http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from alcor.models import Star

logger = logging.getLogger(__name__)


def plot(stars: List[Star],
         *,
         filename: str = 'ugriz.ps',
         figure_size: Tuple[float, float] = (8, 8),
         spacing: float = 0.25,
         ug_label: str = '$u-g$',
         gr_label: str = '$g-r$',
         ri_label: str = '$r-i$',
         iz_label: str = '$i-z$') -> None:
    figure, (subplot_ug_vs_gr,
             subplot_gr_vs_ri,
             subplot_ri_vs_iz) = plt.subplots(nrows=3,
                                              figsize=figure_size)

    # TODO: write a function? relative transformation happens in sampling.py
    # Transformation from UBVRI to ugriz. More info at:
    # Jordi, Grebel & Ammon, 2006, A&A, 460; equations 1-8 and Table 3
    ugriz_ug = [0.75 * (star.u_abs_magnitude - star.b_abs_magnitude) +
                0.77 * (star.b_abs_magnitude - star.v_abs_magnitude) + 0.72
                for star in stars]
    ugriz_gr = [1.646 * (star.v_abs_magnitude - star.r_abs_magnitude) - 0.139
                for star in stars]
    ugriz_ri = [1.007 * (star.r_abs_magnitude - star.i_abs_magnitude) - 0.236
                for star in stars]
    ugriz_iz = [0.577 * (star.r_abs_magnitude - star.i_abs_magnitude) - 0.15
                for star in stars]

    draw_subplot(subplot=subplot_ug_vs_gr,
                 xlabel=gr_label,
                 ylabel=ug_label,
                 x=ugriz_gr,
                 y=ugriz_ug)
    draw_subplot(subplot=subplot_gr_vs_ri,
                 xlabel=ri_label,
                 ylabel=gr_label,
                 x=ugriz_ri,
                 y=ugriz_gr)
    draw_subplot(subplot=subplot_ri_vs_iz,
                 xlabel=iz_label,
                 ylabel=ri_label,
                 x=ugriz_iz,
                 y=ugriz_ri)

    figure.subplots_adjust(hspace=spacing)

    plt.savefig(filename)


def draw_subplot(subplot: Axes,
                 xlabel: str,
                 ylabel: str,
                 x: List[float],
                 y: List[float],
                 color: str = 'b',
                 point_size: float = 0.5,
                 ratio: float = 10 / 13
                 ) -> None:
    subplot.set(xlabel=xlabel,
                ylabel=ylabel)

    subplot.scatter(x=x,
                    y=y,
                    color=color,
                    s=point_size)

    subplot.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

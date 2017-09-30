import logging
import uuid
from typing import (Optional,
                    Tuple,
                    List)

import matplotlib

# More info at
#  http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.colors import Colormap
import numpy as np

from alcor.models import Star
from alcor.services.common import (PECULIAR_SOLAR_VELOCITY_U,
                                   PECULIAR_SOLAR_VELOCITY_V,
                                   PECULIAR_SOLAR_VELOCITY_W)

logger = logging.getLogger(__name__)


def plot(*,
         stars: List[Star],
         axes: str,
         uv_filename: str = 'heatmap_uv.ps',
         uw_filename: str = 'heatmap_uw.ps',
         vw_filename: str = 'heatmap_vw.ps',
         u_label: str = '$U(km/s)$',
         v_label: str = '$V(km/s)$',
         w_label: str = '$W(km/s)$') -> None:
    # TODO: add coordinates
    if axes != 'velocities':
        return

    # TODO: add choosing frame: relative to Sun/LSR. Now it's rel. to LSR
    u_velocities = [float(star.u_velocity) + PECULIAR_SOLAR_VELOCITY_U
                    for star in stars]
    v_velocities = [float(star.v_velocity) + PECULIAR_SOLAR_VELOCITY_V
                    for star in stars]
    w_velocities = [float(star.w_velocity) + PECULIAR_SOLAR_VELOCITY_W
                    for star in stars]

    # TODO: add option of plotting 3 heatmaps in one fig. at the same time
    draw_plot(xlabel=u_label,
              ylabel=v_label,
              xdata=u_velocities,
              ydata=v_velocities,
              filename=uv_filename)
    draw_plot(xlabel=u_label,
              ylabel=w_label,
              xdata=u_velocities,
              ydata=w_velocities,
              filename=uw_filename)
    draw_plot(xlabel=v_label,
              ylabel=w_label,
              xdata=v_velocities,
              ydata=w_velocities,
              filename=vw_filename)


def draw_plot(*,
              xlabel: str,
              ylabel: str,
              xdata: List[float],
              ydata: List[float],
              filename: str,
              figure_size: Tuple[float, float] = (8, 8),
              ratio: float = 10 / 13,
              spacing: float = 0.25,
              figure_grid_height_ratios: List[float] = None,
              bins_count: int = 150,
              vmin: float = 0.01) -> None:
    if figure_grid_height_ratios is None:
        figure_grid_height_ratios = [0.05, 1]

    figure, (colorbar, subplot) = plt.subplots(
            nrows=2,
            figsize=figure_size,
            gridspec_kw={'height_ratios': figure_grid_height_ratios})

    # TODO: add sliders
    subplot.set(xlabel=xlabel,
                ylabel=ylabel)

    heatmap, xedges, yedges = np.histogram2d(x=xdata,
                                             y=ydata,
                                             bins=bins_count)
    extent = [xedges[0], xedges[-1],
              yedges[0], yedges[-1]]

    colormap = colormap_by_name('jet', under='w')

    colorbar_src = subplot.imshow(X=heatmap.T,
                                  cmap=colormap,
                                  vmin=vmin,
                                  extent=extent,
                                  origin='lower')

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    figure.colorbar(mappable=colorbar_src,
                    cax=colorbar,
                    orientation='horizontal')

    subplot.set_aspect(ratio / subplot.get_data_ratio())

    figure.subplots_adjust(hspace=spacing)

    plt.savefig(filename)


def colormap_by_name(name: str,
                     *,
                     under: str = None) -> Colormap:
    colormap = cm.get_cmap(name)
    if under is not None:
        colormap.set_under(under)
    return colormap

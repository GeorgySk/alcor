from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from alcor.models.star import Star
from alcor.services.data_access.reading import fetch


FILENAME = 'heatmap.ps'

FIGURE_SIZE = (8, 8)
DESIRED_DIMENSIONS_RATIO = 10 / 13
SUBPLOTS_SPACING = 0.25

COLORMAP = cm.get_cmap('jet').set_under('w')

U_LABEL = '$U(km/s)$'
V_LABEL = '$V(km/s)$'
W_LABEL = '$W(km/s)$'
VELOCITIES_BINS_COUNT = 75


def plot(*,
         session: Session,
         heatmaps_axes: str) -> None:
    # TODO: Figure out what stars I should fetch (all/last group by time/last N
    # groups by time/selected by ID/marked by some flag(series of simulations))
    stars = fetch_all_stars(session=session)

    figure, (subplot_top,
             subplot_middle,
             subplot_bottom) = plt.subplots(nrows=3,
                                            figsize=FIGURE_SIZE)

    # TODO: add coordinates
    if heatmaps_axes == 'velocities':
        # TODO: add sliders
        subplot_top.set(xlabel=U_LABEL,
                        ylabel=V_LABEL)
        subplot_middle.set(xlabel=U_LABEL,
                           ylabel=W_LABEL)
        subplot_bottom.set(xlabel=V_LABEL,
                           ylabel=W_LABEL)

        velocities_u = [star.velocity_u
                        for star in stars]
        velocities_v = [star.velocity_v
                        for star in stars]
        velocities_w = [star.velocity_w
                        for star in stars]

        heatmap_uv, xedges_uv, yedges_uv = np.histogram2d(
            x=velocities_u,
            y=velocities_v,
            bins=VELOCITIES_BINS_COUNT)
        extent_uv = [xedges_uv[0], xedges_uv[-1],
                     yedges_uv[0], yedges_uv[-1]]
        heatmap_uw, xedges_uw, yedges_uw = np.histogram2d(
            x=velocities_u,
            y=velocities_w,
            bins=VELOCITIES_BINS_COUNT)
        extent_uw = [xedges_uw[0], xedges_uw[-1],
                     yedges_uw[0], yedges_uw[-1]]
        heatmap_vw, xedges_vw, yedges_vw = np.histogram2d(
            x=velocities_v,
            y=velocities_w,
            bins=VELOCITIES_BINS_COUNT)
        extent_vw = [xedges_vw[0], xedges_vw[-1],
                     yedges_vw[0], yedges_vw[-1]]

        subplot_top.imshow(X=heatmap_uv.T,
                           cmap=COLORMAP,
                           vmin=0.01,
                           extent=extent_uv)
        subplot_middle.imshow(X=heatmap_uw.T,
                              cmap=COLORMAP,
                              vmin=0.01,
                              extent=extent_uw)
        subplot_bottom.imshow(X=heatmap_vw.T,
                              cmap=COLORMAP,
                              vmin=0.01,
                              extent=extent_vw)

        # TODO: why does this apply minorticks only to the last subplot?
        plt.minorticks_on()

        subplot_top.xaxis.set_ticks_position('both')
        subplot_top.yaxis.set_ticks_position('both')
        subplot_middle.xaxis.set_ticks_position('both')
        subplot_middle.yaxis.set_ticks_position('both')
        subplot_bottom.xaxis.set_ticks_position('both')
        subplot_bottom.yaxis.set_ticks_position('both')

        subplot_top.set_aspect(DESIRED_DIMENSIONS_RATIO
                               / subplot_top.get_data_ratio())
        subplot_middle.set_aspect(DESIRED_DIMENSIONS_RATIO
                                  / subplot_middle.get_data_ratio())
        subplot_bottom.set_aspect(DESIRED_DIMENSIONS_RATIO
                                  / subplot_bottom.get_data_ratio())

        figure.subplots_adjust(hspace=SUBPLOTS_SPACING)

        plt.savefig(FILENAME)


def fetch_all_stars(*,
                    session: Session):
    query = (Star.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [Star(**record)
            for record in records]

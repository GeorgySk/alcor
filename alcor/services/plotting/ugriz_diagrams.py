import logging

from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.services.data_access.reading import fetch_unprocessed_groups
from alcor.services.processing.service import fetch_stars


logger = logging.getLogger(__name__)


FILENAME = 'ugriz.ps'

FIGURE_SIZE = (8, 8)
DESIRED_DIMENSIONS_RATIO = 10 / 13
SUBPLOTS_SPACING = 0.25

UG_LABEL = '$u-g$'
GR_LABEL = '$g-r$'
RI_LABEL = '$r-i$'
IZ_LABEL = '$i-z$'

POINT_COLOR = 'b'
POINT_SIZE = 0.5


def plot(session: Session) -> None:
    # TODO: Implement other options
    unprocessed_groups = fetch_unprocessed_groups(session=session)
    logger.debug(f'Number of unprocessed groups: {len(unprocessed_groups)}')

    ugriz_ug = []
    ugriz_gr = []
    ugriz_ri = []
    ugriz_iz = []
    for group_count, group in enumerate(unprocessed_groups):
        stars = fetch_stars(group=group,
                            session=session)
        ug = [star.ugriz_ug
              for star in stars]
        ugriz_ug.extend(ug)
        gr = [star.ugriz_gr
              for star in stars]
        ugriz_gr.extend(gr)
        ri = [star.ugriz_ri
              for star in stars]
        ugriz_ri.extend(ri)
        iz = [star.ugriz_iz
              for star in stars]
        ugriz_iz.extend(iz)

    # TODO: do I need to use sharex or sharey attrs?
    figure, (subplot_ug_vs_gr,
             subplot_gr_vs_ri,
             subplot_ri_vs_iz) = plt.subplots(nrows=3,
                                              figsize=FIGURE_SIZE)

    # TODO: find the way to apply limits once for all subplots
    subplot_ug_vs_gr.set(xlabel=GR_LABEL,
                         ylabel=UG_LABEL)
    subplot_gr_vs_ri.set(xlabel=RI_LABEL,
                         ylabel=GR_LABEL)
    subplot_ri_vs_iz.set(xlabel=IZ_LABEL,
                         ylabel=RI_LABEL)
    subplot_ug_vs_gr.scatter(x=ugriz_gr,
                             y=ugriz_ug,
                             color=POINT_COLOR,
                             s=POINT_SIZE)
    subplot_gr_vs_ri.scatter(x=ugriz_ri,
                             y=ugriz_gr,
                             color=POINT_COLOR,
                             s=POINT_SIZE)
    subplot_ri_vs_iz.scatter(x=ugriz_iz,
                             y=ugriz_ri,
                             color=POINT_COLOR,
                             s=POINT_SIZE)

    # TODO: why does this apply minorticks only to the last subplot?
    plt.minorticks_on()

    subplot_ug_vs_gr.xaxis.set_ticks_position('both')
    subplot_ug_vs_gr.yaxis.set_ticks_position('both')
    subplot_gr_vs_ri.xaxis.set_ticks_position('both')
    subplot_gr_vs_ri.yaxis.set_ticks_position('both')
    subplot_ri_vs_iz.xaxis.set_ticks_position('both')
    subplot_ri_vs_iz.yaxis.set_ticks_position('both')

    subplot_ug_vs_gr.set_aspect(DESIRED_DIMENSIONS_RATIO
                                / subplot_ug_vs_gr.get_data_ratio())
    subplot_gr_vs_ri.set_aspect(DESIRED_DIMENSIONS_RATIO
                                / subplot_gr_vs_ri.get_data_ratio())
    subplot_ri_vs_iz.set_aspect(DESIRED_DIMENSIONS_RATIO
                                / subplot_ri_vs_iz.get_data_ratio())

    figure.subplots_adjust(hspace=SUBPLOTS_SPACING)

    plt.savefig(FILENAME)

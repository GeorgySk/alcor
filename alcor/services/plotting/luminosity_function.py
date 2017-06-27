from cassandra.cluster import Session
import matplotlib
# See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for details
# TODO: use this: https://stackoverflow.com/a/37605654/7851470
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alcor.models.luminosity_function import Point
from alcor.services.data_access.reading import fetch


FILENAME = 'luminosity_function.ps'

FIGURE_SIZE = (8, 8)
DESIRED_DIMENSIONS_RATIO = 10 / 13

X_LABEL = '$M_{bol}$'
Y_LABEL = '$\log N (pc^{-3}M_{bol}^{-1})$'

X_LIMITS = [7, 19]
Y_LIMITS = [-6, -2]

LINE_COLOR = 'k'
MARKER = 's'
CAP_SIZE = 5


def plot(*,
         session: Session) -> None:

    # TODO: Implement getting last points by time(ok?)
    graph_points = fetch_all_graph_points(session=session)

    avg_bin_magnitudes = [graph_point.avg_bin_magnitude
                          for graph_point in graph_points]
    stars_count_logarithms = [graph_point.stars_count_logarithm
                              for graph_point in graph_points]
    upper_errorbars = [graph_point.upper_error_bar
                       for graph_point in graph_points]
    lower_errorbars = [graph_point.lower_error_bar
                       for graph_point in graph_points]

    (avg_bin_magnitudes,
     stars_count_logarithms,
     upper_errorbars,
     lower_errorbars) = (_ for _ in zip(*sorted(zip(avg_bin_magnitudes,
                                                    stars_count_logarithms,
                                                    upper_errorbars,
                                                    lower_errorbars))))

    asymmetric_errorbars = [lower_errorbars,
                            upper_errorbars]

    figure, subplot = plt.subplots(figsize=FIGURE_SIZE)

    subplot.set(xlabel=X_LABEL,
                ylabel=Y_LABEL,
                xlim=X_LIMITS,
                ylim=Y_LIMITS)

    subplot.errorbar(x=avg_bin_magnitudes,
                     y=stars_count_logarithms,
                     yerr=asymmetric_errorbars,
                     marker=MARKER,
                     color=LINE_COLOR,
                     capsize=CAP_SIZE)

    plt.minorticks_on()

    subplot.xaxis.set_ticks_position('both')
    subplot.yaxis.set_ticks_position('both')

    subplot.set_aspect(DESIRED_DIMENSIONS_RATIO / subplot.get_data_ratio())

    plt.savefig(FILENAME)


def fetch_all_graph_points(*,
                           session: Session):
    query = (Point.objects.all().limit(None))
    records = fetch(query=query,
                    session=session)
    return [Point(**record)
            for record in records]

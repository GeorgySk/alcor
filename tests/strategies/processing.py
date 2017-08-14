from hypothesis import strategies

from alcor.services.restrictions import FILTRATION_METHODS

filtration_methods = strategies.one_of(map(strategies.just,
                                           FILTRATION_METHODS))

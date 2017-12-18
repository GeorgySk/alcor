import string
from typing import Optional

from hypothesis import strategies
from hypothesis.searchstrategy import SearchStrategy

floats = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-1e15,
                           max_value=1e15)


def positive_floats(max_value: Optional[float] = None) -> SearchStrategy:
    return strategies.floats(min_value=0.,
                             max_value=max_value,
                             allow_nan=False,
                             allow_infinity=False).filter(lambda x: x != 0.)


non_numbers_alphabet = strategies.characters(
        blacklist_characters=string.digits)
non_float_strings = strategies.text(non_numbers_alphabet)

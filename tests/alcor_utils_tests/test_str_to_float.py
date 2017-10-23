from math import isclose

from alcor.utils import str_to_float


def test_str_to_float():
    assert isclose(str_to_float('0.0'), 0.)

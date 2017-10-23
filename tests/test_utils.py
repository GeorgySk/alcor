import math

from alcor.utils import str_to_float


def test_str_to_float(float_value: float,
                      float_value_str: str,
                      non_float_string: str) -> None:
    float_result = str_to_float(float_value_str)
    non_float_result = str_to_float(non_float_string)

    assert isinstance(float_result, float)
    assert isinstance(non_float_result, str)
    assert math.isclose(float_result,
                        float_value)
    assert non_float_result == non_float_string

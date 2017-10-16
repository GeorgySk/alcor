from collections import OrderedDict

from alcor.models import Star


def test_fields_to_copy(star: Star) -> None:
    fields_to_copy = star.fields_to_copy()

    assert isinstance(fields_to_copy, tuple)
    assert all(isinstance(field_name, str) and field_name
               for field_name in fields_to_copy)
    assert fields_to_copy == Star.fields_to_copy()

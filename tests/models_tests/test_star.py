from collections import OrderedDict

from alcor.models import Star


def test_fields_to_copy(star: Star) -> None:
    fields_to_copy = star.fields_to_copy()

    assert isinstance(fields_to_copy, tuple)
    assert all(isinstance(field_name, str) and field_name
               for field_name in fields_to_copy)
    assert fields_to_copy == Star.fields_to_copy()


def test_serialization(star: Star) -> None:
    serialized_star = star.serialize()
    deserialized_star = Star.deserialize(serialized_star)

    assert isinstance(serialized_star, OrderedDict)
    assert all(serialized_star[field_name] == getattr(star, field_name)
               for field_name in Star.fields_to_copy())
    assert all(getattr(deserialized_star, field_name)
               == getattr(star, field_name)
               for field_name in Star.fields_to_copy())

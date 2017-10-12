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


def test_cartesian_coordinates(star: Star) -> None:
    cartesian_coordinates = star.cartesian_coordinates

    assert isinstance(cartesian_coordinates, tuple)
    assert all(isinstance(field_name, float)
               for field_name in cartesian_coordinates)
    assert len(cartesian_coordinates) == 3


def test_max_coordinates_modulus(star: Star) -> None:
    max_coordinates_modulus = star.max_coordinates_modulus

    min_cartesian_coordinates_modulus = min(map(abs,
                                                star.cartesian_coordinates))

    assert isinstance(max_coordinates_modulus, float)
    assert max_coordinates_modulus >= 0.
    assert max_coordinates_modulus >= min_cartesian_coordinates_modulus

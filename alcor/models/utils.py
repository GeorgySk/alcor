from functools import wraps
from typing import (Any,
                    Callable)


def memoize_properties(cls: type) -> type:
    decorator = class_decorator(cls)

    for member_name, member in vars(cls).items():
        if isinstance(member, property):
            original_method = member.fget
            setattr(decorator, member_name, memoized_property(original_method))
    return decorator


def memoized_property(method: Callable[[type], Any]) -> property:
    field_name = '_' + method.__name__

    @wraps(method)
    def memoized(self):
        try:
            return getattr(self, field_name)
        except AttributeError:
            value = method(self)
            setattr(self, field_name, value)
            return value

    return property(memoized)


def class_decorator(cls: type) -> type:
    class decorator(cls):
        pass

    decorator.__doc__ = cls.__doc__
    decorator.__name__ = cls.__name__
    decorator.__qualname__ = cls.__qualname__
    return decorator

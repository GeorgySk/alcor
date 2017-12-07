from hypothesis import strategies


# TODO: is this a good way to make a strategy?
def nonnegative_floats(max_value: float = None):
    return strategies.floats(min_value=0.,
                             max_value=max_value,
                             allow_nan=False,
                             allow_infinity=False)

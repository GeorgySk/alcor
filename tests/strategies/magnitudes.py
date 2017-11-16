from hypothesis import strategies


floats = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-1e15,
                           max_value=1e15)

floats_lists = strategies.lists(elements=floats,
                                min_size=20,
                                max_size=20,
                                unique=True)

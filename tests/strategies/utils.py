import string

from hypothesis import strategies


floats = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-1e15,
                           max_value=1e15)

non_numbers_alphabet = strategies.characters(
        blacklist_characters=string.digits)
non_float_strings = strategies.text(non_numbers_alphabet)

import string

from hypothesis import strategies


floats = strategies.floats(allow_nan=False,
                           allow_infinity=False)

non_numbers_alphabet = strategies.characters(
        blacklist_characters=string.digits)
non_float_strings = strategies.text(non_numbers_alphabet)

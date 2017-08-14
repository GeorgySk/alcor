from collections import Counter

from alcor.models import Star
from alcor.services.group_processing.elimination import check


def test_check(star: Star,
               filtration_method: str) -> None:
    eliminations_counter = Counter()
    empty_eliminations_counter = Counter()

    eliminated = check(star,
                       eliminations_counter=eliminations_counter,
                       filtration_method=filtration_method,
                       min_parallax=float('inf'),
                       min_declination=float('inf'),
                       max_velocity=-float('inf'),
                       min_proper_motion=float('inf'))
    not_eliminated = check(star,
                           eliminations_counter=empty_eliminations_counter,
                           filtration_method='',
                           min_parallax=-float('inf'),
                           min_declination=-float('inf'),
                           max_velocity=float('inf'),
                           min_proper_motion=-float('inf'))

    assert isinstance(eliminated, bool) and isinstance(not_eliminated, bool)
    assert eliminated
    assert eliminations_counter
    assert not not_eliminated
    assert not empty_eliminations_counter

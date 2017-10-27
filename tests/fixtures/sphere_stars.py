from random import random

import pytest


@pytest.fixture(scope='function')
def halo_birth_init_time() -> float:
    return random() * 15.


@pytest.fixture(scope='function')
def halo_stars_formation_time(halo_birth_init_time: float) -> float:
    return random() * (15. - halo_birth_init_time)


@pytest.fixture(scope='function')
def thin_disk_birth_init_time() -> float:
    return random() * 15.


@pytest.fixture(scope='function')
def time_bin() -> float:
    return random()


@pytest.fixture(scope='function')
def time_increment() -> float:
    return random()

from alcor.services.simulations.sphere_stars import (
        halo_star_birth_time,
        thin_disk_star_birth_time,
        thick_disk_star_birth_time)


UNIVERSE_AGE = 14.


def test_halo_star_birth_time(halo_birth_initial_time: float,
                              halo_stars_formation_time: float) -> None:
    birth_time = halo_star_birth_time(
            birth_initial_time=halo_birth_initial_time,
            formation_time=halo_stars_formation_time)

    assert isinstance(birth_time, float)
    assert birth_time >= halo_birth_initial_time
    assert 0. <= birth_time <= birth_time + halo_stars_formation_time


def test_thin_disk_star_birth_time(bin_initial_time: float,
                                   time_increment: float) -> None:
    birth_time = thin_disk_star_birth_time(
            bin_initial_time=bin_initial_time,
            time_increment=time_increment)

    assert isinstance(birth_time, float)
    assert UNIVERSE_AGE >= birth_time >= bin_initial_time


# TODO: find out what to do with RNG
def test_thick_disk_star_birth_time(age: float,
                                    formation_rate_parameter: float,
                                    max_formation_rate: float,
                                    thick_disk_birth_initial_time: float
                                    ) -> None:
    birth_time = thick_disk_star_birth_time(
            age=age,
            formation_rate_parameter=formation_rate_parameter,
            max_formation_rate=max_formation_rate,
            birth_initial_time=thick_disk_birth_initial_time)

    assert isinstance(birth_time, float)
    assert UNIVERSE_AGE >= birth_time >= thick_disk_birth_initial_time

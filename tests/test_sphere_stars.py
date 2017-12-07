from alcor.services.simulations.sphere_stars import (halo_star_birth_time,
                                                     thin_disk_star_birth_time)


def test_halo_star_birth_time(halo_birth_init_time: float,
                              halo_stars_formation_time: float) -> None:
    birth_time = halo_star_birth_time(
            birth_initial_time=halo_birth_init_time,
            formation_time=halo_stars_formation_time)

    assert isinstance(birth_time, float)
    assert birth_time > halo_birth_init_time
    assert birth_time < birth_time + halo_stars_formation_time


def test_thin_disk_star_birth_time(bin_initial_time: float,
                                   time_increment: float) -> None:
    birth_time = thin_disk_star_birth_time(
            bin_initial_time=bin_initial_time,
            time_increment=time_increment)

    assert isinstance(birth_time, float)
    assert birth_time > bin_initial_time

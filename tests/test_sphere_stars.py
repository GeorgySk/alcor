from alcor.services.simulations.sphere_stars import (halo_star_birth_time,
                                                     thin_disk_star_birth_time)


def test_halo_star_birth_time(halo_birth_init_time: float,
                              halo_stars_formation_time: float) -> None:
    birth_time = halo_star_birth_time(
            halo_birth_init_time=halo_birth_init_time,
            halo_stars_formation_time=halo_stars_formation_time)

    assert isinstance(birth_time, float)
    assert birth_time > halo_birth_init_time
    assert birth_time < birth_time + halo_stars_formation_time


def test_thin_disk_star_birth_time(thin_disk_birth_init_time: float,
                                   time_bin: int,
                                   time_increment: float) -> None:
    birth_time = thin_disk_star_birth_time(
            thin_disk_birth_init_time=thin_disk_birth_init_time,
            time_bin=time_bin,
            time_increment=time_increment)

    assert isinstance(birth_time, float)
    assert birth_time > thin_disk_birth_init_time

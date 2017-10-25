import pandas as pd

from alcor.services.simulations.luminosities import (filter_by_max_mass,
                                                     set_metallicities)


def test_filter_by_max_mass(stars_w_progenitor_masses: pd.DataFrame,
                            max_mass: float) -> None:
    filtered_stars = filter_by_max_mass(stars_w_progenitor_masses,
                                        max_mass=max_mass)

    assert isinstance(filtered_stars, pd.DataFrame)
    assert filtered_stars.shape[0] == 1


def test_set_metallicities(stars_w_galactic_disk_types: pd.DataFrame,
                           solar_metallicity: float,
                           subsolar_metallicity: float) -> None:
    set_metallicities(stars_w_galactic_disk_types,
                      subsolar_metallicity=subsolar_metallicity,
                      solar_metallicity=solar_metallicity)
    thin_disk_stars_metallicities = (stars_w_galactic_disk_types[
                stars_w_galactic_disk_types['galactic_disk_type'] == 'thin']
            ['metallicity'])
    thick_disk_stars_metallicities = (stars_w_galactic_disk_types[
                                         stars_w_galactic_disk_types[
                                             'galactic_disk_type'] == 'thick']
                                      ['metallicity'])
    halo_disk_stars_metallicities = (stars_w_galactic_disk_types[
                                         stars_w_galactic_disk_types[
                                             'galactic_disk_type'] == 'halo']
                                     ['metallicity'])

    assert isinstance(stars_w_galactic_disk_types, pd.DataFrame)
    assert stars_w_galactic_disk_types.shape[0] > 0
    assert 'galactic_disk_type' in stars_w_galactic_disk_types.columns
    assert 'metallicity' in stars_w_galactic_disk_types.columns
    assert (thin_disk_stars_metallicities == solar_metallicity).all()
    assert (thick_disk_stars_metallicities == solar_metallicity).all()
    assert (halo_disk_stars_metallicities == subsolar_metallicity).all()

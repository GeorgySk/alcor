import pandas as pd

from alcor.services.simulations.luminosities import filter_by_max_mass


def test_filter_by_max_mass(stars: pd.DataFrame,
                            max_mass: float) -> None:
    filtered_stars = filter_by_max_mass(stars,
                                        max_mass=max_mass)

    assert isinstance(filtered_stars, pd.DataFrame)
    assert filtered_stars.shape[0] == 1

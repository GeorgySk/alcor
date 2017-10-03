import logging
from typing import (Callable,
                    List)

from alcor.models import (Group,
                          Star)
from alcor.types import (BolometricIndexType,
                         StarBolometricIndexType,
                         StarsBinsType)

logger = logging.getLogger(__name__)

# Motion of the Sun taken with respect to the Local Standard of Rest (LSR)
# according to (Schönrich R., Binney J., Dehnen W., 2010, MNRAS, 403, 1829)
PECULIAR_SOLAR_VELOCITY_U = -11
PECULIAR_SOLAR_VELOCITY_V = 12
PECULIAR_SOLAR_VELOCITY_W = 7

FILTRATION_METHODS = ['raw', 'full', 'restricted']

STARS_SPECTRAL_TYPES = {0: 'DA',
                        1: 'DB',
                        2: 'ONe'}


def group_output_file_name(group: Group,
                           *,
                           extension='.res',
                           file_name_length: int = 5) -> str:
    group_id_str = str(group.id)
    base_name = group_id_str[:file_name_length]
    return ''.join([base_name, extension])


def bolometric_indexer(*,
                       min_magnitude: float,
                       stars_bin_size: float) -> BolometricIndexType:
    def bolometric_index(magnitude: float) -> int:
        magnitude_amplitude = magnitude - min_magnitude
        return int(magnitude_amplitude / stars_bin_size)

    return bolometric_index


def star_bolometric_indexer(bolometric_index: Callable[[float], int]
                            ) -> StarBolometricIndexType:
    def star_bolometric_index(star: Star) -> int:
        return bolometric_index(star.bolometric_magnitude)

    return star_bolometric_index


def stars_packer(*,
                 stars_bins_count: int,
                 star_bolometric_index: StarBolometricIndexType
                 ) -> Callable[[List[Star]],
                               StarsBinsType]:
    def pack_stars(stars: List[Star]) -> StarsBinsType:
        res = [[] for _ in range(stars_bins_count)]
        for star in stars:
            index = star_bolometric_index(star)
            if stars_bins_count > index >= 0:
                res[index].append(star)
            else:
                logger.warning(
                        f'Bolometric magnitude {star.bolometric_magnitude} '
                        'is out of bounds '
                        f'for star with id {star.id}.')
        return res

    return pack_stars

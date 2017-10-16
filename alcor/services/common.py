import logging

from alcor.models import Group

logger = logging.getLogger(__name__)

# Motion of the Sun taken with respect to the Local Standard of Rest (LSR)
# according to (Schönrich R., Binney J., Dehnen W., 2010, MNRAS, 403, 1829)
PECULIAR_SOLAR_VELOCITY_U = -11
PECULIAR_SOLAR_VELOCITY_V = 12
PECULIAR_SOLAR_VELOCITY_W = 7

FILTRATION_METHODS = ['raw', 'full', 'restricted']


def group_output_file_name(group: Group,
                           *,
                           extension='.res',
                           file_name_length: int = 5) -> str:
    group_id_str = str(group.id)
    base_name = group_id_str[:file_name_length]
    return ''.join([base_name, extension])

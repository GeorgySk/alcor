import logging

logger = logging.getLogger(__name__)

# Motion of the Sun taken with respect to the Local Standard of Rest (LSR)
# according to (Schönrich R., Binney J., Dehnen W., 2010, MNRAS, 403, 1829)
PECULIAR_SOLAR_VELOCITY_U = -11
PECULIAR_SOLAR_VELOCITY_V = 12
PECULIAR_SOLAR_VELOCITY_W = 7

FILTRATION_METHODS = ['raw', 'full', 'restricted']

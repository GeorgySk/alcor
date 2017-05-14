import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model

STAR_PARAMETERS_NAMES = ['luminosity',
                         'proper_motion',
                         'declination_galactocentric_distance',
                         'gr_photometry',
                         'rz_photometry',
                         'v_photometry',
                         'velocity_u',
                         'velocity_v',
                         'velocity_w']


class Star(Model):
    __table_name__ = 'stars'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    luminosity = Decimal(required=True)
    proper_motion = Decimal(required=True)
    declination_galactocentric_distance = Decimal(required=True)
    gr_photometry = Decimal(required=True)
    rz_photometry = Decimal(required=True)
    v_photometry = Decimal(required=True)
    velocity_u = Decimal(required=True)
    velocity_v = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)

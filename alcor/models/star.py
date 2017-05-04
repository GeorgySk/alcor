import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Float,
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
    luminosity = Float(required=True)
    proper_motion = Float(required=True)
    declination_galactocentric_distance = Float(required=True)
    gr_photometry = Float(required=True)
    rz_photometry = Float(required=True)
    v_photometry = Float(required=True)
    velocity_u = Float(required=True)
    velocity_v = Float(required=True)
    velocity_w = Float(required=True)
    updated_timestamp = DateTime(default=datetime.now)

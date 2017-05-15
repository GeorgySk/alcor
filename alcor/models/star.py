import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model

STAR_PARAMETERS_NAMES = ['luminosity',
                         'proper_motion',
                         'declination',
                         'galactocentric_distance',
                         'go_photometry'
                         'gr_photometry',
                         'rz_photometry',
                         'v_photometry',
                         'velocity_u',
                         'velocity_v',
                         'velocity_w']


class Star(Model):
    def __init__(self,
                 parameters: list):
        self.luminosity = parameters[0]
        self.proper_motion = parameters[1]
        self.declination = parameters[2]
        self.galactocentric_distance = parameters[3]
        self.go_photometry = parameters[4]
        self.gr_photometry = parameters[5]
        self.rz_photometry = parameters[6]
        self.v_photometry = parameters[7]
        self.velocity_u = parameters[8]
        self.velocity_v = parameters[9]
        self.velocity_w = parameters[10]

    __table_name__ = 'stars'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    luminosity = Decimal(required=True)
    proper_motion = Decimal(required=True)
    declination = Decimal(required=True)
    galactocentric_distance = Decimal(required=True)
    go_photometry = Decimal(required=True)
    gr_photometry = Decimal(required=True)
    rz_photometry = Decimal(required=True)
    v_photometry = Decimal(required=True)
    velocity_u = Decimal(required=True)
    velocity_v = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)

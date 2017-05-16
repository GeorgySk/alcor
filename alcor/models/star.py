import uuid
from datetime import datetime
from math import (cos,
                  sin)

from astropy import units as u
from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime,
                                         Integer)
from cassandra.cqlengine.models import Model

ASTRONOMICAL_UNIT = 4.74 * u.km / u.s

STAR_PARAMETERS_NAMES = ['luminosity',
                         'proper_motion',
                         'proper_motion_component_b',
                         'proper_motion_component_l',
                         'proper_motion_component_vr',
                         'right_ascension',
                         'declination',
                         'galactocentric_distance',
                         'galactocentric_coordinate_b',
                         'galactocentric_coordinate_l',
                         'go_photometry',
                         'gr_photometry',
                         'rz_photometry',
                         'v_photometry',
                         'velocity_u',
                         'velocity_v',
                         'velocity_w',
                         'spectral_type']


class Star(Model):
    __table_name__ = 'stars'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    luminosity = Decimal(required=True)
    proper_motion = Decimal(required=True)
    proper_motion_component_b = Decimal(required=True)
    proper_motion_component_l = Decimal(required=True)
    proper_motion_component_vr = Decimal(required=True)
    right_ascension = Decimal(required=True)
    declination = Decimal(required=True)
    galactocentric_distance = Decimal(required=True)
    galactocentric_coordinate_b = Decimal(required=True)
    galactocentric_coordinate_l = Decimal(required=True)
    go_photometry = Decimal(required=True)
    gr_photometry = Decimal(required=True)
    rz_photometry = Decimal(required=True)
    v_photometry = Decimal(required=True)
    velocity_u = Decimal(required=True)
    velocity_v = Decimal(required=True)
    velocity_w = Decimal(required=True)
    spectral_type = Integer(required=True)
    updated_timestamp = DateTime(default=datetime.now)

    def set_radial_velocity_to_zero(self) -> None:
        distance_in_pc = self.galactocentric_distance * 10e3

        a1 = (-ASTRONOMICAL_UNIT
              * cos(self.galactocentric_coordinate_b)
              * sin(self.galactocentric_coordinate_l))
        b1 = (-ASTRONOMICAL_UNIT
              * sin(self.galactocentric_coordinate_b)
              * cos(self.galactocentric_coordinate_l))
        self.velocity_u = (a1 * self.proper_motion_component_l * distance_in_pc
                           + b1 * self.proper_motion_component_b
                              * distance_in_pc)

        a2 = (ASTRONOMICAL_UNIT
              * cos(self.galactocentric_coordinate_b)
              * cos(self.galactocentric_coordinate_l))
        b2 = (-ASTRONOMICAL_UNIT
              * sin(self.galactocentric_coordinate_b)
              * sin(self.galactocentric_coordinate_l))
        self.velocity_v = (a2 * self.proper_motion_component_l * distance_in_pc
                           + b2 * self.proper_motion_component_b
                              * distance_in_pc)

        b3 = ASTRONOMICAL_UNIT * cos(
            self.galactocentric_coordinate_b)
        self.velocity_w = (b3 * self.proper_motion_component_b
                           * distance_in_pc)

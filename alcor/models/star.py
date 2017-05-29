import decimal
import uuid
from datetime import datetime
from math import (cos,
                  sin,
                  pi,
                  asin,
                  atan)
from typing import Tuple
from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime,
                                         Integer)
from cassandra.cqlengine.models import Model

ASTRONOMICAL_UNIT = 4.74
DEC_GPOLE = 27.128336 * pi / 180.0
RA_GPOLE = 192.859508 * pi / 180.0
AUX_ANGLE = 122.932 * pi / 180.0

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

    @property
    def bolometric_magnitude(self) -> float:
        # TODO: find out the meaning of the following constants
        return decimal.Decimal(2.5) * self.luminosity + decimal.Decimal(4.75)

    @property
    def coordinate_x(self) -> float:
        return float(self.to_cartesian_from_equatorial()[0])

    @property
    def coordinate_y(self) -> float:
        return float(self.to_cartesian_from_equatorial()[1])

    @property
    def coordinate_z(self) -> float:
        return float(self.to_cartesian_from_equatorial()[2])

    def to_cartesian_from_equatorial(self) -> Tuple[float,
                                                    float,
                                                    float]:
        ra = float(self.right_ascension)
        dec = float(self.declination)
        distance = float(self.galactocentric_distance)

        b = (asin(cos(dec) * cos(DEC_GPOLE) * cos(ra - RA_GPOLE)
                  + sin(dec) * sin(DEC_GPOLE)))
        x = sin(dec) - sin(b) * sin(DEC_GPOLE)
        y = cos(dec) * sin(ra - RA_GPOLE) * cos(DEC_GPOLE)
        l = atan(x / y) + AUX_ANGLE - pi / 2.0
        if y < 0. and x > 0.:
            l += pi
        elif x < 0. and y < 0.:
            l += pi
        elif x < 0. and y > 0.:
            l += 2.0 * pi
        elif l > 2.0 * pi:
            l -= 2.0 * pi
        coordinate_x = distance * cos(b) * cos(l)
        coordinate_y = distance * cos(b) * sin(l)
        coordinate_z = distance * sin(b)
        return (coordinate_x,
                coordinate_y,
                coordinate_z)

    def set_radial_velocity_to_zero(self) -> None:
        # TODO: implement pc/kpc units
        galactocentric_distance = float(self.galactocentric_distance)
        galactocentric_coordinate_b = float(self.galactocentric_coordinate_b)
        galactocentric_coordinate_l = float(self.galactocentric_coordinate_l)
        proper_motion_component_b = float(self.proper_motion_component_b)
        proper_motion_component_l = float(self.proper_motion_component_l)

        distance_in_pc = galactocentric_distance * 1e3

        a1 = (-ASTRONOMICAL_UNIT * cos(galactocentric_coordinate_b)
              * sin(galactocentric_coordinate_l))
        b1 = (-ASTRONOMICAL_UNIT * sin(galactocentric_coordinate_b)
              * cos(galactocentric_coordinate_l))
        self.velocity_u = ((a1 * proper_motion_component_l
                           + b1 * proper_motion_component_b)
                           * distance_in_pc)

        a2 = (ASTRONOMICAL_UNIT * cos(galactocentric_coordinate_b)
              * cos(galactocentric_coordinate_l))
        b2 = (-ASTRONOMICAL_UNIT * sin(galactocentric_coordinate_b)
              * sin(galactocentric_coordinate_l))
        self.velocity_v = ((a2 * proper_motion_component_l
                           + b2 * proper_motion_component_b)
                           * distance_in_pc)

        b3 = ASTRONOMICAL_UNIT * cos(galactocentric_coordinate_b)
        self.velocity_w = (b3 * proper_motion_component_b
                           * distance_in_pc)

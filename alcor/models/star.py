import decimal
import uuid
from datetime import datetime
from math import (cos,
                  sin)
from typing import Tuple

from astropy import units as u
from astropy.coordinates.sky_coordinate import SkyCoord
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

    @property
    def bolometric_magnitude(self) -> float:
        # TODO: find out the meaning of the following constants
        return decimal.Decimal(2.5) * self.luminosity + decimal.Decimal(4.75)

    @property
    def coordinate_x(self) -> float:
        try:
            return self.coordinate_x
        except AttributeError:
            return self.to_cartesian_from_equatorial[0]

    @property
    def coordinate_y(self) -> float:
        try:
            return self.coordinate_y
        except AttributeError:
            return self.to_cartesian_from_equatorial[1]

    @property
    def coordinate_z(self) -> float:
        try:
            return self.coordinate_z
        except AttributeError:
            return self.to_cartesian_from_equatorial[2]

    def to_cartesian_from_equatorial(self) -> Tuple[Decimal,
                                                    Decimal,
                                                    Decimal]:
        equatorial_coordinates = SkyCoord(
            ra=self.right_ascension * u.degree,
            dec=self.declination * u.degree,
            distance=self.distance * u.kpc)
        return (equatorial_coordinates.cartesian.x,
                equatorial_coordinates.cartesian.y,
                equatorial_coordinates.cartesian.z)

    def set_radial_velocity_to_zero(self) -> None:
        # TODO: implement pc/kpc units
        distance_in_pc = self.galactocentric_distance * 10e3

        a1 = (-ASTRONOMICAL_UNIT * cos(self.galactocentric_coordinate_b)
              * sin(self.galactocentric_coordinate_l))
        b1 = (-ASTRONOMICAL_UNIT * sin(self.galactocentric_coordinate_b)
              * cos(self.galactocentric_coordinate_l))
        self.velocity_u = ((a1 * self.proper_motion_component_l
                           + b1 * self.proper_motion_component_b)
                           * distance_in_pc)

        a2 = (ASTRONOMICAL_UNIT * cos(self.galactocentric_coordinate_b)
              * cos(self.galactocentric_coordinate_l))
        b2 = (-ASTRONOMICAL_UNIT * sin(self.galactocentric_coordinate_b)
              * sin(self.galactocentric_coordinate_l))
        self.velocity_v = ((a2 * self.proper_motion_component_l
                           + b2 * self.proper_motion_component_b)
                           * distance_in_pc)

        b3 = ASTRONOMICAL_UNIT * cos(self.galactocentric_coordinate_b)
        self.velocity_w = (b3 * self.proper_motion_component_b
                           * distance_in_pc)

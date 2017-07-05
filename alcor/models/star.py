import uuid
from decimal import Decimal
from math import (cos,
                  sin,
                  pi,
                  asin,
                  atan)
from typing import Tuple

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (BigInteger,
                                     Integer,
                                     Float,
                                     DateTime)

from .base import Base

ASTRONOMICAL_UNIT = 4.74
DEC_GPOLE = 27.128336 * pi / 180.
RA_GPOLE = 192.859508 * pi / 180.
AUX_ANGLE = 122.932 * pi / 180.

STAR_PARAMETERS_NAMES = ['luminosity',
                         'proper_motion',
                         'proper_motion_component_b',
                         'proper_motion_component_l',
                         'proper_motion_component_vr',
                         'right_ascension',
                         'declination',
                         'galactic_distance',
                         'galactic_latitude',
                         'galactic_longitude',
                         'go_photometry',
                         'gr_photometry',
                         'rz_photometry',
                         'v_photometry',
                         'velocity_u',
                         'velocity_v',
                         'velocity_w',
                         'spectral_type',
                         'disk_belonging']


class Star(Base):
    __tablename__ = 'stars'

    id = Column(BigInteger(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    luminosity = Column(Float(asdecimal=True),
                        nullable=False)
    proper_motion = Column(Float(asdecimal=True),
                           nullable=False)
    proper_motion_component_b = Column(Float(asdecimal=True),
                                       nullable=False)
    proper_motion_component_l = Column(Float(asdecimal=True),
                                       nullable=False)
    proper_motion_component_vr = Column(Float(asdecimal=True),
                                        nullable=False)
    right_ascension = Column(Float(asdecimal=True),
                             nullable=False)
    declination = Column(Float(asdecimal=True),
                         nullable=False)
    galactic_distance = Column(Float(asdecimal=True),
                               nullable=False)
    galactic_coordinate_b = Column(Float(asdecimal=True),
                                   nullable=False)
    galactic_coordinate_l = Column(Float(asdecimal=True),
                                   nullable=False)
    go_photometry = Column(Float(asdecimal=True),
                           nullable=False)
    gr_photometry = Column(Float(asdecimal=True),
                           nullable=False)
    rz_photometry = Column(Float(asdecimal=True),
                           nullable=False)
    v_photometry = Column(Float(asdecimal=True),
                          nullable=False)
    velocity_u = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_v = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_w = Column(Float(asdecimal=True),
                        nullable=False)
    spectral_type = Column(Integer(),
                           nullable=False)
    disk_belonging = Column(Integer(),
                            nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 luminosity: Decimal,
                 proper_motion: Decimal,
                 proper_motion_component_b: Decimal,
                 proper_motion_component_l: Decimal,
                 proper_motion_component_vr: Decimal,
                 right_ascension: Decimal,
                 declination: Decimal,
                 galactic_distance: Decimal,
                 galactic_coordinate_b: Decimal,
                 galactic_coordinate_l: Decimal,
                 go_photometry: Decimal,
                 gr_photometry: Decimal,
                 rz_photometry: Decimal,
                 v_photometry: Decimal,
                 velocity_u: Decimal,
                 velocity_v: Decimal,
                 velocity_w: Decimal,
                 spectral_type: int,
                 disk_belonging: int):
        self.group_id = group_id
        self.luminosity = luminosity
        self.proper_motion = proper_motion
        self.proper_motion_component_b = proper_motion_component_b
        self.proper_motion_component_l = proper_motion_component_l
        self.proper_motion_component_vr = proper_motion_component_vr
        self.right_ascension = right_ascension
        self.declination = declination
        self.galactic_distance = galactic_distance
        self.galactic_coordinate_b = galactic_coordinate_b
        self.galactic_coordinate_l = galactic_coordinate_l
        self.go_photometry = go_photometry
        self.gr_photometry = gr_photometry
        self.rz_photometry = rz_photometry
        self.v_photometry = v_photometry
        self.velocity_u = velocity_u
        self.velocity_v = velocity_v
        self.velocity_w = velocity_w
        self.spectral_type = spectral_type
        self.disk_belonging = disk_belonging


    @property
    def bolometric_magnitude(self) -> float:
        # TODO: find out the meaning of the following constants
        return Decimal(2.5) * self.luminosity + Decimal(4.75)

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
        right_ascension = float(self.right_ascension)
        declination = float(self.declination)
        # TODO: should it be centered at the Sun?
        distance = float(self.galactic_distance)

        latitude = (asin(cos(declination) * cos(DEC_GPOLE)
                         * cos(right_ascension - RA_GPOLE)
                         + sin(declination) * sin(DEC_GPOLE)))
        x = sin(declination) - sin(latitude) * sin(DEC_GPOLE)
        y = cos(declination) * sin(right_ascension - RA_GPOLE) * cos(DEC_GPOLE)
        longitude = atan(x / y) + AUX_ANGLE - pi / 2.
        if x > 0. and 0. > y or x <= 0. and y <= 0.:
            longitude += pi
        coordinate_x = distance * cos(latitude) * cos(longitude)
        coordinate_y = distance * cos(latitude) * sin(longitude)
        coordinate_z = distance * sin(latitude)
        return (coordinate_x,
                coordinate_y,
                coordinate_z)

    def set_radial_velocity_to_zero(self) -> None:
        # TODO: implement pc/kpc units
        galactic_distance = float(self.galactic_distance)
        galactic_latitude = float(self.galactic_latitude)
        galactic_longitude = float(self.galactic_longitude)
        proper_motion_component_b = float(self.proper_motion_component_b)
        proper_motion_component_l = float(self.proper_motion_component_l)

        distance_in_pc = galactic_distance * 1e3

        a1 = (-ASTRONOMICAL_UNIT * cos(galactic_latitude)
              * sin(galactic_longitude))
        b1 = (-ASTRONOMICAL_UNIT * sin(galactic_latitude)
              * cos(galactic_longitude))
        self.velocity_u = ((a1 * proper_motion_component_l
                            + b1 * proper_motion_component_b)
                           * distance_in_pc)

        a2 = (ASTRONOMICAL_UNIT * cos(galactic_latitude)
              * cos(galactic_longitude))
        b2 = (-ASTRONOMICAL_UNIT * sin(galactic_latitude)
              * sin(galactic_longitude))
        self.velocity_v = ((a2 * proper_motion_component_l
                            + b2 * proper_motion_component_b)
                           * distance_in_pc)

        b3 = ASTRONOMICAL_UNIT * cos(galactic_latitude)
        self.velocity_w = (b3 * proper_motion_component_b
                           * distance_in_pc)

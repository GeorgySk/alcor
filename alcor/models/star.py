import enum
import uuid
# TODO: figure out if to use Decimal at all or smth else from scipy package
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
                                     DateTime, Enum)

from .base import Base

ASTRONOMICAL_UNIT = 4.74
DEC_GPOLE = 27.128336 * pi / 180.
RA_GPOLE = 192.859508 * pi / 180.
AUX_ANGLE = 122.932 * pi / 180.

SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE = 4.75

# Description of ugriz color scheme:
# https://en.wikipedia.org/wiki/Photometric_system
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
                         'ugriz_g_apparent',
                         'ugriz_ug',
                         'ugriz_gr',
                         'ugriz_ri',
                         'ugriz_iz',
                         'v_photometry',
                         'velocity_u',
                         'velocity_v',
                         'velocity_w',
                         'spectral_type',
                         'disk_belonging']


class GalacticDiskEnum(enum.Enum):
    thin = 1
    thick = 2


class Star(Base):
    __tablename__ = 'stars'

    id = Column(BigInteger(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    luminosity = Column(Float(asdecimal=True),
                        nullable=True)
    proper_motion = Column(Float(asdecimal=True),
                           nullable=True)
    proper_motion_component_b = Column(Float(asdecimal=True),
                                       nullable=True)
    proper_motion_component_l = Column(Float(asdecimal=True),
                                       nullable=True)
    proper_motion_component_vr = Column(Float(asdecimal=True),
                                        nullable=True)
    right_ascension = Column(Float(asdecimal=True),
                             nullable=True)
    declination = Column(Float(asdecimal=True),
                         nullable=True)
    galactic_distance = Column(Float(asdecimal=True),
                               nullable=True)
    galactic_latitude = Column(Float(asdecimal=True),
                               nullable=True)
    galactic_longitude = Column(Float(asdecimal=True),
                                nullable=True)
    ugriz_g_apparent = Column(Float(asdecimal=True),
                              nullable=True)
    ugriz_ug = Column(Float(asdecimal=True),
                      nullable=True)
    ugriz_gr = Column(Float(asdecimal=True),
                      nullable=True)
    ugriz_ri = Column(Float(asdecimal=True),
                      nullable=True)
    ugriz_iz = Column(Float(asdecimal=True),
                      nullable=True)
    v_photometry = Column(Float(asdecimal=True),
                          nullable=True)
    velocity_u = Column(Float(asdecimal=True),
                        nullable=True)
    velocity_v = Column(Float(asdecimal=True),
                        nullable=True)
    velocity_w = Column(Float(asdecimal=True),
                        nullable=True)
    spectral_type = Column(Integer(),
                           nullable=True)
    disk_belonging = Column(Enum(GalacticDiskEnum),
                            nullable=True)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 luminosity: float = None,
                 proper_motion: float = None,
                 proper_motion_component_b: float = None,
                 proper_motion_component_l: float = None,
                 proper_motion_component_vr: float = None,
                 right_ascension: float = None,
                 declination: float = None,
                 galactic_distance: float = None,
                 galactic_latitude: float = None,
                 galactic_longitude: float = None,
                 ugriz_g_apparent: float = None,
                 ugriz_ug: float = None,
                 ugriz_gr: float = None,
                 ugriz_ri: float = None,
                 ugriz_iz: float = None,
                 v_photometry: float = None,
                 velocity_u: float = None,
                 velocity_v: float = None,
                 velocity_w: float = None,
                 spectral_type: int = None,
                 disk_belonging: int = None):
        self.group_id = group_id
        self.luminosity = luminosity
        self.proper_motion = proper_motion
        self.proper_motion_component_b = proper_motion_component_b
        self.proper_motion_component_l = proper_motion_component_l
        self.proper_motion_component_vr = proper_motion_component_vr
        self.right_ascension = right_ascension
        self.declination = declination
        self.galactic_distance = galactic_distance
        self.galactic_latitude = galactic_latitude
        self.galactic_longitude = galactic_longitude
        self.ugriz_g_apparent = ugriz_g_apparent
        self.ugriz_ug = ugriz_ug
        self.ugriz_gr = ugriz_gr
        self.ugriz_ri = ugriz_ri
        self.ugriz_iz = ugriz_iz
        self.v_photometry = v_photometry
        self.velocity_u = velocity_u
        self.velocity_v = velocity_v
        self.velocity_w = velocity_w
        self.spectral_type = spectral_type
        self.disk_belonging = disk_belonging

    @property
    def bolometric_magnitude(self) -> float:
        # More info at
        # https://en.wikipedia.org/wiki/Absolute_magnitude#Bolometric_magnitude
        return (2.5 * float(self.luminosity)
                + SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE)

    @property
    def coordinate_x(self) -> float:
        return float(self.to_cartesian_from_equatorial()[0])

    @property
    def coordinate_y(self) -> float:
        return float(self.to_cartesian_from_equatorial()[1])

    @property
    def coordinate_z(self) -> float:
        return float(self.to_cartesian_from_equatorial()[2])

    @property
    def ugriz_rz(self) -> float:
        return float(self.ugriz_ri) + float(self.ugriz_iz)

    def to_cartesian_from_equatorial(self) -> Tuple[float,
                                                    float,
                                                    float]:
        right_ascension = float(self.right_ascension)
        declination = float(self.declination)
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

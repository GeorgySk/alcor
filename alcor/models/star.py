import enum
import inspect
import uuid
from collections import OrderedDict
from math import (radians,
                  cos,
                  sin,
                  pi,
                  asin,
                  atan)
from typing import (Any,
                    Dict,
                    Tuple)

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (BigInteger,
                                     Integer,
                                     Float,
                                     DateTime,
                                     Enum)

from .base import Base
from .utils import memoize_properties

ASTRONOMICAL_UNIT = 4.74
DEC_GPOLE = radians(27.128336)
RA_GPOLE = radians(192.859508)
AUX_ANGLE = radians(122.932)

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
                         'distance',
                         'galactic_latitude',
                         'galactic_longitude',
                         'ugriz_g_apparent',
                         'ugriz_ug',
                         'ugriz_gr',
                         'ugriz_ri',
                         'ugriz_iz',
                         'v_photometry',
                         'u_velocity',
                         'v_velocity',
                         'w_velocity',
                         'spectral_type',
                         'galactic_disk_type']


class GalacticDiskType(enum.IntEnum):
    thin = 1
    thick = 2


@memoize_properties
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
    distance = Column(Float(asdecimal=True),
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
    u_velocity = Column(Float(asdecimal=True),
                        nullable=True)
    v_velocity = Column(Float(asdecimal=True),
                        nullable=True)
    w_velocity = Column(Float(asdecimal=True),
                        nullable=True)
    # TODO: make it Enum, DA - 0, DB - 1, ONe - 2
    spectral_type = Column(Integer(),
                           nullable=True)
    galactic_disk_type = Column(Enum(GalacticDiskType),
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
                 distance: float = None,
                 galactic_latitude: float = None,
                 galactic_longitude: float = None,
                 ugriz_g_apparent: float = None,
                 ugriz_ug: float = None,
                 ugriz_gr: float = None,
                 ugriz_ri: float = None,
                 ugriz_iz: float = None,
                 v_photometry: float = None,
                 u_velocity: float = None,
                 v_velocity: float = None,
                 w_velocity: float = None,
                 spectral_type: int = None,
                 galactic_disk_type: GalacticDiskType = None):
        self.id = None
        self.group_id = group_id
        self.luminosity = luminosity
        self.proper_motion = proper_motion
        self.proper_motion_component_b = proper_motion_component_b
        self.proper_motion_component_l = proper_motion_component_l
        self.proper_motion_component_vr = proper_motion_component_vr
        self.right_ascension = right_ascension
        self.declination = declination
        self.distance = distance
        self.galactic_latitude = galactic_latitude
        self.galactic_longitude = galactic_longitude
        self.ugriz_g_apparent = ugriz_g_apparent
        self.ugriz_ug = ugriz_ug
        self.ugriz_gr = ugriz_gr
        self.ugriz_ri = ugriz_ri
        self.ugriz_iz = ugriz_iz
        self.v_photometry = v_photometry
        self.u_velocity = u_velocity
        self.v_velocity = v_velocity
        self.w_velocity = w_velocity
        self.spectral_type = spectral_type
        self.galactic_disk_type = galactic_disk_type

    @property
    def bolometric_magnitude(self) -> float:
        # More info at
        # https://en.wikipedia.org/wiki/Absolute_magnitude#Bolometric_magnitude
        return (2.5 * float(self.luminosity)
                + SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE)

    @property
    def x_coordinate(self) -> float:
        return float(self.cartesian_coordinates[0])

    @property
    def y_coordinate(self) -> float:
        return float(self.cartesian_coordinates[1])

    @property
    def z_coordinate(self) -> float:
        return float(self.cartesian_coordinates[2])

    @property
    def max_coordinates_modulus(self) -> float:
        return max(abs(self.x_coordinate),
                   abs(self.y_coordinate),
                   abs(self.z_coordinate))

    @property
    def ugriz_rz(self) -> float:
        return float(self.ugriz_ri) + float(self.ugriz_iz)

    @property
    def cartesian_coordinates(self) -> Tuple[float, float, float]:
        right_ascension = float(self.right_ascension)
        declination = float(self.declination)
        distance = float(self.distance)

        latitude = (asin(cos(declination) * cos(DEC_GPOLE)
                         * cos(right_ascension - RA_GPOLE)
                         + sin(declination) * sin(DEC_GPOLE)))
        x = sin(declination) - sin(latitude) * sin(DEC_GPOLE)
        y = cos(declination) * sin(right_ascension - RA_GPOLE) * cos(DEC_GPOLE)
        longitude = atan(x / y) + AUX_ANGLE - pi / 2.
        if x > 0. and 0. > y or x <= 0. and y <= 0.:
            longitude += pi

        x_coordinate = distance * cos(latitude) * cos(longitude)
        y_coordinate = distance * cos(latitude) * sin(longitude)
        z_coordinate = distance * sin(latitude)
        return x_coordinate, y_coordinate, z_coordinate

    def modify(self, **fields: Any) -> 'Star':
        serialized_star = self.serialize()
        serialized_star.update(fields)
        return self.deserialize(serialized_star)

    def serialize(self) -> Dict[str, Any]:
        return OrderedDict((field_name, getattr(self, field_name))
                           for field_name in self.fields_to_copy())

    @classmethod
    def deserialize(cls, serialized: Dict[str, Any]) -> 'Star':
        serialized = dict(serialized)
        star_id = serialized.pop('id')
        star = cls(**serialized)
        star.id = star_id
        return star

    # TODO: memoize this
    @classmethod
    def fields_to_copy(cls) -> Tuple[str, ...]:
        initializer_signature = inspect.signature(cls.__init__)
        parameters = dict(initializer_signature.parameters)
        parameters.pop('self')
        return ('id',) + tuple(parameters)


def set_radial_velocity_to_zero(star: Star) -> Star:
    # TODO: implement pc/kpc units
    distance = float(star.distance)
    galactic_latitude = float(star.galactic_latitude)
    galactic_longitude = float(star.galactic_longitude)
    proper_motion_component_b = float(star.proper_motion_component_b)
    proper_motion_component_l = float(star.proper_motion_component_l)

    distance_in_pc = distance * 1e3

    a1 = (-ASTRONOMICAL_UNIT * cos(galactic_latitude)
          * sin(galactic_longitude))
    b1 = (-ASTRONOMICAL_UNIT * sin(galactic_latitude)
          * cos(galactic_longitude))
    u_velocity = ((a1 * proper_motion_component_l
                   + b1 * proper_motion_component_b)
                  * distance_in_pc)

    a2 = (ASTRONOMICAL_UNIT * cos(galactic_latitude)
          * cos(galactic_longitude))
    b2 = (-ASTRONOMICAL_UNIT * sin(galactic_latitude)
          * sin(galactic_longitude))
    v_velocity = ((a2 * proper_motion_component_l
                   + b2 * proper_motion_component_b)
                  * distance_in_pc)

    b3 = ASTRONOMICAL_UNIT * cos(galactic_latitude)
    w_velocity = (b3 * proper_motion_component_b
                  * distance_in_pc)
    return star.modify(u_velocity=u_velocity,
                       v_velocity=v_velocity,
                       w_velocity=w_velocity)

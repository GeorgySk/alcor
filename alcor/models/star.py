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
STAR_PARAMETERS_NAMES = {'mass',
                         'luminosity',
                         'proper_motion',
                         'proper_motion_component_b',
                         'proper_motion_component_l',
                         'proper_motion_component_vr',
                         'right_ascension',
                         'declination',
                         'r_galactocentric',
                         'th_galactocentric',
                         'z_coordinate',
                         'right_ascension_proper_motion',
                         'declination_proper_motion',
                         'distance',
                         'galactic_latitude',
                         'galactic_longitude',
                         'j_abs_magnitude',
                         'b_abs_magnitude',
                         'v_abs_magnitude',
                         'r_abs_magnitude',
                         'i_abs_magnitude',
                         'u_velocity',
                         'v_velocity',
                         'w_velocity',
                         'birth_time',
                         'spectral_type',
                         'galactic_disk_type'}


class GalacticDiskType(enum.IntEnum):
    thin = 1
    thick = 2
    halo = 3


@memoize_properties
class Star(Base):
    __tablename__ = 'stars'

    id = Column(BigInteger(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    mass = Column(Float(),
                  nullable=True)
    luminosity = Column(Float(),
                        nullable=True)
    r_galactocentric = Column(Float(),
                              nullable=True)
    th_galactocentric = Column(Float(),
                               nullable=True)
    z_coordinate = Column(Float(),
                          nullable=True)
    proper_motion = Column(Float(),
                           nullable=True)
    proper_motion_component_b = Column(Float(),
                                       nullable=True)
    proper_motion_component_l = Column(Float(),
                                       nullable=True)
    proper_motion_component_vr = Column(Float(),
                                        nullable=True)
    right_ascension = Column(Float(),
                             nullable=True)
    declination = Column(Float(),
                         nullable=True)
    right_ascension_proper_motion = Column(Float(),
                                           nullable=True)
    declination_proper_motion = Column(Float(),
                                       nullable=True)
    distance = Column(Float(),
                      nullable=True)
    galactic_latitude = Column(Float(),
                               nullable=True)
    galactic_longitude = Column(Float(),
                                nullable=True)
    j_abs_magnitude = Column(Float(),
                             nullable=True)
    b_abs_magnitude = Column(Float(),
                             nullable=True)
    r_abs_magnitude = Column(Float(),
                             nullable=True)
    v_abs_magnitude = Column(Float(),
                             nullable=True)
    i_abs_magnitude = Column(Float(),
                             nullable=True)
    u_velocity = Column(Float(),
                        nullable=True)
    v_velocity = Column(Float(),
                        nullable=True)
    w_velocity = Column(Float(),
                        nullable=True)
    birth_time = Column(Float(),
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
                 mass: float = None,
                 luminosity: float = None,
                 proper_motion: float = None,
                 r_galactocentric: float = None,
                 th_galactocentric: float = None,
                 z_coordinate: float = None,
                 proper_motion_component_b: float = None,
                 proper_motion_component_l: float = None,
                 proper_motion_component_vr: float = None,
                 right_ascension: float = None,
                 declination: float = None,
                 right_ascension_proper_motion: float = None,
                 declination_proper_motion: float = None,
                 distance: float = None,
                 galactic_latitude: float = None,
                 galactic_longitude: float = None,
                 j_abs_magnitude: float = None,
                 b_abs_magnitude: float = None,
                 r_abs_magnitude: float = None,
                 v_abs_magnitude: float = None,
                 i_abs_magnitude: float = None,
                 u_velocity: float = None,
                 v_velocity: float = None,
                 w_velocity: float = None,
                 birth_time: float = None,
                 spectral_type: int = None,
                 galactic_disk_type: GalacticDiskType = None):
        self.id = None
        self.group_id = group_id
        self.mass = mass
        self.luminosity = luminosity
        self.proper_motion = proper_motion
        self.r_galactocentric = r_galactocentric
        self.th_galactocentric = th_galactocentric
        self.z_coordinate = z_coordinate
        self.proper_motion_component_b = proper_motion_component_b
        self.proper_motion_component_l = proper_motion_component_l
        self.proper_motion_component_vr = proper_motion_component_vr
        self.right_ascension = right_ascension
        self.declination = declination
        self.right_ascension_proper_motion = right_ascension_proper_motion
        self.declination_proper_motion = declination_proper_motion
        self.distance = distance
        self.galactic_latitude = galactic_latitude
        self.galactic_longitude = galactic_longitude
        self.j_abs_magnitude = j_abs_magnitude
        self.b_abs_magnitude = b_abs_magnitude
        self.r_abs_magnitude = r_abs_magnitude
        self.v_abs_magnitude = v_abs_magnitude
        self.i_abs_magnitude = i_abs_magnitude
        self.u_velocity = u_velocity
        self.v_velocity = v_velocity
        self.w_velocity = w_velocity
        self.birth_time = birth_time
        self.spectral_type = spectral_type
        self.galactic_disk_type = galactic_disk_type

    @property
    def bolometric_magnitude(self) -> float:
        # More info at
        # https://en.wikipedia.org/wiki/Absolute_magnitude#Bolometric_magnitude
        return 2.5 * self.luminosity + SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE

    @property
    def x_coordinate(self) -> float:
        return self.cartesian_coordinates[0]

    @property
    def y_coordinate(self) -> float:
        return self.cartesian_coordinates[1]

    @property
    def max_coordinates_modulus(self) -> float:
        return max(abs(self.x_coordinate),
                   abs(self.y_coordinate),
                   abs(self.z_coordinate))

    @property
    def cartesian_coordinates(self) -> Tuple[float, float, float]:
        right_ascension = self.right_ascension
        declination = self.declination
        distance = self.distance

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
        fields_to_copy = cls.fields_to_copy()
        serialized = {key: value
                      for key, value in serialized.items()
                      if key in fields_to_copy}
        star_id = serialized.pop('id', None)
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

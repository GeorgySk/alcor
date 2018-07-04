
import enum
import uuid

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (BigInteger,
                                     Float,
                                     DateTime,
                                     Enum)

from .base import Base

# Description of ugriz color scheme:
# https://en.wikipedia.org/wiki/Photometric_system
STAR_PARAMETERS_NAMES = {'mass',
                         'progenitor_mass',
                         'luminosity',
                         'proper_motion',
                         'proper_motion_in_latitude',
                         'proper_motion_in_longitude',
                         'radial_velocity',
                         'right_ascension',
                         'declination',
                         'r_cylindrical',
                         'theta_cylindrical',
                         'z_coordinate',
                         'right_ascension_proper_motion',
                         'declination_proper_motion',
                         'distance',
                         'galactic_latitude',
                         'galactic_longitude',
                         'color_u',
                         'color_b',
                         'color_v',
                         'color_r',
                         'color_i',
                         'color_j',
                         'cooling_time',
                         'effective_temperature',
                         'metallicity',
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


class SpectralType(enum.IntEnum):
    DA = 0
    DB = 1
    ONe = 2


class Star(Base):
    __tablename__ = 'stars'

    id = Column(BigInteger(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    mass = Column(Float(),
                  nullable=True)
    progenitor_mass = Column(Float(),
                             nullable=True)
    luminosity = Column(Float(),
                        nullable=True)
    r_cylindrical = Column(Float(),
                           nullable=True)
    theta_cylindrical = Column(Float(),
                               nullable=True)
    z_coordinate = Column(Float(),
                          nullable=True)
    proper_motion = Column(Float(),
                           nullable=True)
    proper_motion_in_latitude = Column(Float(),
                                       nullable=True)
    proper_motion_in_longitude = Column(Float(),
                                        nullable=True)
    radial_velocity = Column(Float(),
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
    color_u = Column(Float(),
                     nullable=True)
    color_b = Column(Float(),
                     nullable=True)
    color_v = Column(Float(),
                     nullable=True)
    color_r = Column(Float(),
                     nullable=True)
    color_i = Column(Float(),
                     nullable=True)
    color_j = Column(Float(),
                     nullable=True)
    cooling_time = Column(Float(),
                          nullable=True)
    effective_temperature = Column(Float(),
                                   nullable=True)
    metallicity = Column(Float(),
                         nullable=True)
    u_velocity = Column(Float(),
                        nullable=True)
    v_velocity = Column(Float(),
                        nullable=True)
    w_velocity = Column(Float(),
                        nullable=True)
    birth_time = Column(Float(),
                        nullable=True)
    # TODO: fix these enums
    # spectral_type = Column(Enum(SpectralType),
    #                        nullable=True)
    # galactic_disk_type = Column(Enum(GalacticDiskType),
    #                             nullable=True)
    spectral_type = Column(BigInteger(),
                           nullable=True)
    galactic_disk_type = Column(BigInteger(),
                                nullable=True)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 mass: float = None,
                 progenitor_mass: float = None,
                 luminosity: float = None,
                 proper_motion: float = None,
                 r_cylindrical: float = None,
                 theta_cylindrical: float = None,
                 z_coordinate: float = None,
                 proper_motion_in_latitude: float = None,
                 proper_motion_in_longitude: float = None,
                 radial_velocity: float = None,
                 right_ascension: float = None,
                 declination: float = None,
                 right_ascension_proper_motion: float = None,
                 declination_proper_motion: float = None,
                 distance: float = None,
                 galactic_latitude: float = None,
                 galactic_longitude: float = None,
                 color_u: float = None,
                 color_b: float = None,
                 color_v: float = None,
                 color_r: float = None,
                 color_i: float = None,
                 color_j: float = None,
                 cooling_time: float = None,
                 effective_temperature: float = None,
                 metallicity: float = None,
                 u_velocity: float = None,
                 v_velocity: float = None,
                 w_velocity: float = None,
                 birth_time: float = None,
                 # TODO: fix these enums
                 # spectral_type: SpectralType = None,
                 # galactic_disk_type: GalacticDiskType = None):
                 spectral_type: int = None,
                 galactic_disk_type: int = None):
        self.id = None
        self.group_id = group_id
        self.mass = mass
        self.progenitor_mass = progenitor_mass
        self.luminosity = luminosity
        self.proper_motion = proper_motion
        self.r_cylindrical = r_cylindrical
        self.theta_cylindrical = theta_cylindrical
        self.z_coordinate = z_coordinate
        self.proper_motion_in_latitude = proper_motion_in_latitude
        self.proper_motion_in_longitude = proper_motion_in_longitude
        self.radial_velocity = radial_velocity
        self.right_ascension = right_ascension
        self.declination = declination
        self.right_ascension_proper_motion = right_ascension_proper_motion
        self.declination_proper_motion = declination_proper_motion
        self.distance = distance
        self.galactic_latitude = galactic_latitude
        self.galactic_longitude = galactic_longitude
        self.color_u = color_u
        self.color_b = color_b
        self.color_v = color_v
        self.color_r = color_r
        self.color_i = color_i
        self.color_j = color_j
        self.cooling_time = cooling_time
        self.effective_temperature = effective_temperature
        self.metallicity = metallicity
        self.u_velocity = u_velocity
        self.v_velocity = v_velocity
        self.w_velocity = w_velocity
        self.birth_time = birth_time
        self.spectral_type = spectral_type
        self.galactic_disk_type = galactic_disk_type

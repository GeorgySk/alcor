import uuid

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (Integer,
                                     DateTime)

from .base import Base


class StarsCounter(Base):
    __tablename__ = 'stars_eliminations_counters'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    raw = Column(Integer(),
                 nullable=False)
    by_parallax = Column(Integer(),
                         nullable=False)
    by_declination = Column(Integer(),
                            nullable=False)
    by_velocity = Column(Integer(),
                         nullable=False)
    by_proper_motion = Column(Integer(),
                              nullable=False)
    by_reduced_proper_motion = Column(Integer(),
                                      nullable=False)
    by_apparent_magnitude = Column(Integer(),
                                   nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 raw: int,
                 by_parallax: int = 0,
                 by_declination: int = 0,
                 by_velocity: int = 0,
                 by_proper_motion: int = 0,
                 by_reduced_proper_motion: int = 0,
                 by_apparent_magnitude: int = 0):
        self.group_id = group_id
        self.raw = raw
        self.by_parallax = by_parallax
        self.by_declination = by_declination
        self.by_velocity = by_velocity
        self.by_proper_motion = by_proper_motion
        self.by_reduced_proper_motion = by_reduced_proper_motion
        self.by_apparent_magnitude = by_apparent_magnitude

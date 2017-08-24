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
    parallax = Column(Integer(),
                      nullable=False)
    declination = Column(Integer(),
                         nullable=False)
    velocity = Column(Integer(),
                      nullable=False)
    proper_motion = Column(Integer(),
                           nullable=False)
    reduced_proper_motion = Column(Integer(),
                                   nullable=False)
    apparent_magnitude = Column(Integer(),
                                nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 raw: int,
                 parallax: int,
                 declination: int,
                 velocity: int,
                 proper_motion: int,
                 reduced_proper_motion: int,
                 apparent_magnitude: int):
        self.group_id = group_id
        self.raw = raw
        self.parallax = parallax
        self.declination = declination
        self.velocity = velocity
        self.proper_motion = proper_motion
        self.reduced_proper_motion = reduced_proper_motion
        self.apparent_magnitude = apparent_magnitude

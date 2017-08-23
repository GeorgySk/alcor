import uuid

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (Integer,
                                     Float,
                                     DateTime)

from alcor.models.base import Base


class Bin(Base):
    __tablename__ = 'velocities_vs_magnitudes_bins'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    avg_magnitude = Column(Float(asdecimal=True),
                           nullable=False)
    avg_u_velocity = Column(Float(asdecimal=True),
                            nullable=False)
    avg_v_velocity = Column(Float(asdecimal=True),
                            nullable=False)
    avg_w_velocity = Column(Float(asdecimal=True),
                            nullable=False)
    u_velocity_std = Column(Float(asdecimal=True),
                            nullable=False)
    v_velocity_std = Column(Float(asdecimal=True),
                            nullable=False)
    w_velocity_std = Column(Float(asdecimal=True),
                            nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_u_velocity: float,
                 avg_v_velocity: float,
                 avg_w_velocity: float,
                 u_velocity_std: float,
                 v_velocity_std: float,
                 w_velocity_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_u_velocity = avg_u_velocity
        self.avg_v_velocity = avg_v_velocity
        self.avg_w_velocity = avg_w_velocity
        self.u_velocity_std = u_velocity_std
        self.v_velocity_std = v_velocity_std
        self.w_velocity_std = w_velocity_std


class LepineCaseUBin(Base):
    __tablename__ = 'lepine_case_u_velocities_vs_magnitudes_bins'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    avg_magnitude = Column(Float(asdecimal=True),
                           nullable=False)
    avg_u_velocity = Column(Float(asdecimal=True),
                            nullable=False)
    u_velocity_std = Column(Float(asdecimal=True),
                            nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_u_velocity: float,
                 u_velocity_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_u_velocity = avg_u_velocity
        self.u_velocity_std = u_velocity_std


class LepineCaseVBin(Base):
    __tablename__ = 'lepine_case_v_velocities_vs_magnitudes_bins'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    avg_magnitude = Column(Float(asdecimal=True),
                           nullable=False)
    avg_v_velocity = Column(Float(asdecimal=True),
                            nullable=False)
    v_velocity_std = Column(Float(asdecimal=True),
                            nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_v_velocity: float,
                 v_velocity_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_v_velocity = avg_v_velocity
        self.v_velocity_std = v_velocity_std


class LepineCaseWBin(Base):
    __tablename__ = 'lepine_case_w_velocities_vs_magnitudes_bins'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    avg_magnitude = Column(Float(asdecimal=True),
                           nullable=False)
    avg_w_velocity = Column(Float(asdecimal=True),
                            nullable=False)
    w_velocity_std = Column(Float(asdecimal=True),
                            nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_w_velocity: float,
                 w_velocity_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_w_velocity = avg_w_velocity
        self.w_velocity_std = w_velocity_std

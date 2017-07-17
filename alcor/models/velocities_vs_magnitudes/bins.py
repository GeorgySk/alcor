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
    avg_velocity_u = Column(Float(asdecimal=True),
                            nullable=False)
    avg_velocity_v = Column(Float(asdecimal=True),
                            nullable=False)
    avg_velocity_w = Column(Float(asdecimal=True),
                            nullable=False)
    velocity_u_std = Column(Float(asdecimal=True),
                            nullable=False)
    velocity_v_std = Column(Float(asdecimal=True),
                            nullable=False)
    velocity_w_std = Column(Float(asdecimal=True),
                            nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_velocity_u: float,
                 avg_velocity_v: float,
                 avg_velocity_w: float,
                 velocity_u_std: float,
                 velocity_v_std: float,
                 velocity_w_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_velocity_u = avg_velocity_u
        self.avg_velocity_v = avg_velocity_v
        self.avg_velocity_w = avg_velocity_w
        self.velocity_u_std = velocity_u_std
        self.velocity_v_std = velocity_v_std
        self.velocity_w_std = velocity_w_std


class LepineCaseUBin(Base):
    __tablename__ = 'lepine_case_u_velocities_vs_magnitudes_bins'

    id = Column(Integer(), primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    avg_magnitude = Column(Float(asdecimal=True), nullable=False)
    avg_velocity_u = Column(Float(asdecimal=True), nullable=False)
    velocity_u_std = Column(Float(asdecimal=True), nullable=False)
    updated_timestamp = Column(DateTime(), server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_velocity_u: float,
                 velocity_u_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_velocity_u = avg_velocity_u
        self.velocity_u_std = velocity_u_std


class LepineCaseVBin(Base):
    __tablename__ = 'lepine_case_v_velocities_vs_magnitudes_bins'

    id = Column(Integer(), primary_key=True)
    group_id = Column(UUID(as_uuid=True), nullable=False)
    avg_magnitude = Column(Float(asdecimal=True), nullable=False)
    avg_velocity_v = Column(Float(asdecimal=True), nullable=False)
    velocity_v_std = Column(Float(asdecimal=True), nullable=False)
    updated_timestamp = Column(DateTime(), server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_velocity_v: float,
                 velocity_v_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_velocity_v = avg_velocity_v
        self.velocity_v_std = velocity_v_std


class LepineCaseWBin(Base):
    __tablename__ = 'lepine_case_w_velocities_vs_magnitudes_bins'

    id = Column(Integer(), primary_key=True)
    group_id = Column(UUID(as_uuid=True), nullable=False)
    avg_magnitude = Column(Float(asdecimal=True), nullable=False)
    avg_velocity_w = Column(Float(asdecimal=True), nullable=False)
    velocity_w_std = Column(Float(asdecimal=True), nullable=False)
    updated_timestamp = Column(DateTime(), server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_magnitude: float,
                 avg_velocity_w: float,
                 velocity_w_std: float):
        self.group_id = group_id
        self.avg_magnitude = avg_magnitude
        self.avg_velocity_w = avg_velocity_w
        self.velocity_w_std = velocity_w_std

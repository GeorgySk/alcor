import uuid

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (Integer,
                                     Float,
                                     DateTime)

from alcor.models.base import Base


class Cloud(Base):
    __tablename__ = 'velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(asdecimal=True),
                                  nullable=False)
    velocity_u = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_v = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_w = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 bolometric_magnitude: float,
                 velocity_u: float,
                 velocity_v: float,
                 velocity_w: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.velocity_u = velocity_u
        self.velocity_v = velocity_v
        self.velocity_w = velocity_w


class LepineCaseUCloud(Base):
    __tablename__ = 'lepine_case_u_velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(asdecimal=True),
                                  nullable=False)
    velocity_u = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 bolometric_magnitude: float,
                 velocity_u: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.velocity_u = velocity_u


class LepineCaseVCloud(Base):
    __tablename__ = 'lepine_case_v_velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(asdecimal=True),
                                  nullable=False)
    velocity_v = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 bolometric_magnitude: float,
                 velocity_v: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.velocity_v = velocity_v


class LepineCaseWCloud(Base):
    __tablename__ = 'lepine_case_w_velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(asdecimal=True),
                                  nullable=False)
    velocity_w = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 bolometric_magnitude: float,
                 velocity_w: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.velocity_w = velocity_w

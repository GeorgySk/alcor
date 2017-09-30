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
    bolometric_magnitude = Column(Float(),
                                  nullable=False)
    u_velocity = Column(Float(),
                        nullable=False)
    v_velocity = Column(Float(),
                        nullable=False)
    w_velocity = Column(Float(),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 bolometric_magnitude: float,
                 u_velocity: float,
                 v_velocity: float,
                 w_velocity: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.u_velocity = u_velocity
        self.v_velocity = v_velocity
        self.w_velocity = w_velocity


class LepineCaseUCloud(Base):
    __tablename__ = 'lepine_case_u_velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(),
                                  nullable=False)
    u_velocity = Column(Float(),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 *,
                 group_id: uuid.UUID = None,
                 bolometric_magnitude: float,
                 u_velocity: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.u_velocity = u_velocity


class LepineCaseVCloud(Base):
    __tablename__ = 'lepine_case_v_velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(),
                                  nullable=False)
    v_velocity = Column(Float(),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 *,
                 group_id: uuid.UUID = None,
                 bolometric_magnitude: float,
                 v_velocity: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.v_velocity = v_velocity


class LepineCaseWCloud(Base):
    __tablename__ = 'lepine_case_w_velocities_vs_magnitudes_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    bolometric_magnitude = Column(Float(),
                                  nullable=False)
    w_velocity = Column(Float(),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 *,
                 group_id: uuid.UUID = None,
                 bolometric_magnitude: float,
                 w_velocity: float):
        self.group_id = group_id
        self.bolometric_magnitude = bolometric_magnitude
        self.w_velocity = w_velocity

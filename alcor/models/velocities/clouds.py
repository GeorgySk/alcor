import uuid
from decimal import Decimal

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (Integer,
                                     Float,
                                     DateTime)

from alcor.models.base import Base


class Cloud(Base):
    __tablename__ = 'velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True))
    u_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    v_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    w_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 u_velocity: Decimal,
                 v_velocity: Decimal,
                 w_velocity: Decimal):
        self.group_id = group_id
        self.u_velocity = u_velocity
        self.v_velocity = v_velocity
        self.w_velocity = w_velocity


class LepineCaseUVCloud(Base):
    __tablename__ = 'lepine_case_uv_velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    u_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    v_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 u_velocity: Decimal,
                 v_velocity: Decimal):
        self.group_id = group_id
        self.u_velocity = u_velocity
        self.v_velocity = v_velocity


class LepineCaseUWCloud(Base):
    __tablename__ = 'lepine_case_uw_velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    u_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    w_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 u_velocity: Decimal,
                 w_velocity: Decimal):
        self.group_id = group_id
        self.u_velocity = u_velocity
        self.w_velocity = w_velocity


class LepineCaseVWCloud(Base):
    __tablename__ = 'lepine_case_vw_velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    v_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    w_velocity = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 v_velocity: Decimal,
                 w_velocity: Decimal):
        self.group_id = group_id
        self.v_velocity = v_velocity
        self.w_velocity = w_velocity

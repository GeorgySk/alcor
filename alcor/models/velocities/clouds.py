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
                 velocity_u: Decimal,
                 velocity_v: Decimal,
                 velocity_w: Decimal):
        self.group_id = group_id
        self.velocity_u = velocity_u
        self.velocity_v = velocity_v
        self.velocity_w = velocity_w


class LepineCaseUVCloud(Base):
    __tablename__ = 'lepine_case_uv_velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    velocity_u = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_v = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 velocity_u: Decimal,
                 velocity_v: Decimal):
        self.group_id = group_id
        self.velocity_u = velocity_u
        self.velocity_v = velocity_v


class LepineCaseUWCloud(Base):
    __tablename__ = 'lepine_case_uw_velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    velocity_u = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_w = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 velocity_u: Decimal,
                 velocity_w: Decimal):
        self.group_id = group_id
        self.velocity_u = velocity_u
        self.velocity_w = velocity_w


class LepineCaseVWCloud(Base):
    __tablename__ = 'lepine_case_vw_velocities_clouds'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    velocity_v = Column(Float(asdecimal=True),
                        nullable=False)
    velocity_w = Column(Float(asdecimal=True),
                        nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 velocity_v: Decimal,
                 velocity_w: Decimal):
        self.group_id = group_id
        self.velocity_v = velocity_v
        self.velocity_w = velocity_w

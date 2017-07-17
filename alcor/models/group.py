import uuid

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import (Boolean,
                                     DateTime)

from .base import Base


class Group(Base):
    __tablename__ = 'groups'

    id = Column(UUID(as_uuid=True),
                primary_key=True)
    processed = Column(Boolean(),
                       default=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self, id: uuid.UUID):
        self.id = id

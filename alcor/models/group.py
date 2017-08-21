import uuid
from typing import Optional

from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.sqltypes import DateTime

from .base import Base


class Group(Base):
    __tablename__ = 'groups'
    id = Column(UUID(as_uuid=True),
                primary_key=True)
    original_id = Column(UUID(as_uuid=True),
                         default=None)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 id: uuid.UUID,
                 original_id: Optional[uuid.UUID]):
        self.id = id
        self.original_id = original_id

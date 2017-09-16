from sqlalchemy.schema import Column
from sqlalchemy.sql.schema import ForeignKey
from sqlalchemy.sql.sqltypes import BigInteger

from .base import Base
from .star import Star


class ProcessedStars(Base):
    __tablename__ = 'processed_stars_associations'
    original_id = Column(BigInteger(),
                         ForeignKey(Star.id),
                         primary_key=True)
    processed_id = Column(BigInteger(),
                          ForeignKey(Star.id),
                          primary_key=True)

    def __init__(self,
                 original_id: int,
                 processed_id: int):
        self.original_id = original_id
        self.processed_id = processed_id

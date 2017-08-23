from sqlalchemy.schema import Column
from sqlalchemy.sql.schema import ForeignKey
from sqlalchemy.sql.sqltypes import BigInteger

from .base import Base
from .star import Star


class ProcessedStarsAssociation(Base):
    __tablename__ = 'processed_stars_associations'
    original_star_id = Column(BigInteger(),
                              ForeignKey(Star.id),
                              primary_key=True)
    processed_star_id = Column(BigInteger(),
                               ForeignKey(Star.id),
                               primary_key=True)

    def __init__(self,
                 original_star_id: int,
                 processed_star_id: int):
        self.original_star_id = original_star_id
        self.processed_star_id = processed_star_id

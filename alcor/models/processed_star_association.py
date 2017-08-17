from sqlalchemy.schema import Column
from sqlalchemy.sql.sqltypes import BigInteger

from .base import Base


# TODO: fix this
class ProcessedStarAssociation(Base):
    __tablename__ = 'processed_stars_associations'
    id = Column(BigInteger(),
                primary_key=True)
    original_star_id = Column(BigInteger(),
                              default=None)

    def __init__(self,
                 original_star_id: int):
        self.original_star_id = original_star_id

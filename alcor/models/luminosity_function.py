import uuid

from sqlalchemy import func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.schema import Column
from sqlalchemy.sql.sqltypes import (Integer,
                                     Float,
                                     DateTime)

from .base import Base


class Point(Base):
    __tablename__ = 'luminosity_function_graph'

    id = Column(Integer(),
                primary_key=True)
    group_id = Column(UUID(as_uuid=True),
                      nullable=False)
    avg_bin_magnitude = Column(Float(asdecimal=True),
                               nullable=False)
    stars_count_logarithm = Column(Float(asdecimal=True),
                                   nullable=False)
    upper_error_bar = Column(Float(asdecimal=True),
                             nullable=False)
    lower_error_bar = Column(Float(asdecimal=True),
                             nullable=False)
    updated_timestamp = Column(DateTime(),
                               server_default=func.now())

    def __init__(self,
                 group_id: uuid.UUID,
                 avg_bin_magnitude: float,
                 stars_count_logarithm: float,
                 upper_error_bar: float,
                 lower_error_bar: float):
        self.group_id = group_id
        self.avg_bin_magnitude = avg_bin_magnitude
        self.stars_count_logarithm = stars_count_logarithm
        self.upper_error_bar = upper_error_bar
        self.lower_error_bar = lower_error_bar

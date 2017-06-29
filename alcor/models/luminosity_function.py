import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model


class Point(Model):
    __table_name__ = 'luminosity_function_graph'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    avg_bin_magnitude = Decimal(required=True)
    stars_count_logarithm = Decimal(required=True)
    upper_error_bar = Decimal(required=True)
    lower_error_bar = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now,
                                 primary_key=True,
                                 clustering_order='DESC')

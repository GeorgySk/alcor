import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Text,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model


class CParameter(Model):
    __table_name__ = 'simulations_parameters'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True,
                    index=True)
    name = Text(required=True)
    value = Decimal(required=True)
    created_timestamp = DateTime(default=datetime.now)

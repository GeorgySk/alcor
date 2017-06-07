import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Boolean,
                                         DateTime)
from cassandra.cqlengine.models import Model


class Group(Model):
    __table_name__ = 'groups'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    processed = Boolean(index=True,
                        default=False)
    updated_timestamp = DateTime(default=datetime.now)

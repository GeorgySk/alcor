import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Integer,
                                         DateTime)
from cassandra.cqlengine.models import Model


class CStarsCounter(Model):
    __table_name__ = 'stars_eliminations_counters'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    raw = Integer(required=True)
    by_parallax = Integer(required=True)
    by_declination = Integer(required=True)
    by_velocity = Integer(required=True)
    by_proper_motion = Integer(required=True)
    by_reduced_proper_motion = Integer(required=True)
    by_apparent_magnitude = Integer(required=True)
    updated_timestamp = DateTime(default=datetime.now)

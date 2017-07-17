import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model


class CCloud(Model):
    __table_name__ = 'velocities_vs_magnitudes_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    bolometric_magnitude = Decimal(required=True)
    velocity_u = Decimal(required=True)
    velocity_v = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class CLepineCaseUCloud(Model):
    __table_name__ = 'lepine_case_u_velocities_vs_magnitudes_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    bolometric_magnitude = Decimal(required=True)
    velocity_u = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class CLepineCaseVCloud(Model):
    __table_name__ = 'lepine_case_v_velocities_vs_magnitudes_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    bolometric_magnitude = Decimal(required=True)
    velocity_v = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class CLepineCaseWCloud(Model):
    __table_name__ = 'lepine_case_w_velocities_vs_magnitudes_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    bolometric_magnitude = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)

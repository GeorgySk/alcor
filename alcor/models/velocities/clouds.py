import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model


class Cloud(Model):
    __table_name__ = 'velocities_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    velocity_u = Decimal(required=True)
    velocity_v = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class LepineCaseUVCloud(Model):
    __table_name__ = 'lepine_case_uv_velocities_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    velocity_u = Decimal(required=True)
    velocity_v = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class LepineCaseUWCloud(Model):
    __table_name__ = 'lepine_case_uw_velocities_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    velocity_u = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class LepineCaseVWCloud(Model):
    __table_name__ = 'lepine_case_vw_velocities_clouds'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    velocity_v = Decimal(required=True)
    velocity_w = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)

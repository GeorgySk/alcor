import uuid
from datetime import datetime

from cassandra.cqlengine.columns import (UUID,
                                         Decimal,
                                         DateTime)
from cassandra.cqlengine.models import Model


class Bin(Model):
    __table_name__ = 'velocities_vs_magnitudes_bins'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    avg_magnitude = Decimal(required=True)
    avg_velocity_u = Decimal(required=True)
    avg_velocity_v = Decimal(required=True)
    avg_velocity_w = Decimal(required=True)
    velocity_u_std = Decimal(required=True)
    velocity_v_std = Decimal(required=True)
    velocity_w_std = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class LepineCaseUBin(Model):
    __table_name__ = 'lepine_case_u_velocities_vs_magnitudes_bins'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    avg_magnitude = Decimal(required=True)
    avg_velocity_u = Decimal(required=True)
    velocity_u_std = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class LepineCaseVBin(Model):
    __table_name__ = 'lepine_case_v_velocities_vs_magnitudes_bins'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    avg_magnitude = Decimal(required=True)
    avg_velocity_v = Decimal(required=True)
    velocity_v_std = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)


class LepineCaseWBin(Model):
    __table_name__ = 'lepine_case_w_velocities_vs_magnitudes_bins'

    id = UUID(primary_key=True,
              default=uuid.uuid4)
    group_id = UUID(required=True)
    avg_magnitude = Decimal(required=True)
    avg_velocity_w = Decimal(required=True)
    velocity_w_std = Decimal(required=True)
    updated_timestamp = DateTime(default=datetime.now)

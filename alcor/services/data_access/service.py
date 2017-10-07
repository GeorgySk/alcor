import uuid
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import Group
from alcor.models.base import Base


def fetch_all(model: Base,
              *,
              session: Session) -> List[Base]:
    query = session.query(model)
    return query.all()


def fetch_last_groups(*,
                      limit: int,
                      session: Session) -> List[Group]:
    query = (session.query(Group)
             .order_by(Group.updated_timestamp.desc())
             .limit(limit))
    return query.all()


def fetch_group_by_id(group_id: uuid.UUID,
                      *,
                      session: Session) -> Group:
    query = (session.query(Group)
             .filter(Group.id == group_id))
    return query.one()

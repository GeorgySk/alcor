import uuid
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models.base import Base
from alcor.models.group import Group
from alcor.models.star import Star


def fetch_all(model: Base,
              *,
              session: Session) -> List[Base]:
    query = session.query(model)
    return query.all()


def fetch_group_stars(*,
                      group_id: uuid.UUID,
                      session: Session) -> List[Star]:
    query = (session.query(Star)
             .filter(Star.group_id == group_id))
    return query.all()


def fetch_unprocessed_groups(session: Session) -> List[Group]:
    query = (session.query(Group)
             .filter(Group.original_id.is_(None)))
    return query.all()


def fetch_last_groups(*,
                      limit: int,
                      session: Session) -> List[Group]:
    query = (session.query(Group)
             .order_by(Group.updated_timestamp.desc())
             .limit(limit))
    return query.all()


def fetch_group_by_id(*,
                      group_id: uuid.UUID,
                      session: Session) -> Group:
    query = (session.query(Group)
             .filter(Group.id == group_id))
    return query.one()

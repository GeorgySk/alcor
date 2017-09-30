import uuid
from typing import Optional

from sqlalchemy.orm.session import Session

from alcor.services import stars_group
from alcor.services.data_access import (fetch_unprocessed_groups,
                                        fetch_last_groups,
                                        fetch_group_by_id)


def run(*,
        last_groups_count: Optional[int],
        unprocessed_groups: bool,
        group_id: Optional[uuid.UUID],
        session: Session) -> None:
    if unprocessed_groups:
        groups = fetch_unprocessed_groups(session=session)
    elif last_groups_count:
        groups = fetch_last_groups(limit=last_groups_count,
                                   session=session)
    else:
        groups = [fetch_group_by_id(group_id=group_id,
                                    session=session)]
    for group in groups:
        stars_group.process(
                group=group,
                session=session)

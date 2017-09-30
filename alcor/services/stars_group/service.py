import uuid
from typing import (Iterator,
                    List)

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star,
                          associations)


def process(*,
            group: Group,
            session: Session) -> None:

    original_id = group.id
    processed_group_id = uuid.uuid4()
    processed_group = Group(id=processed_group_id,
                            original_id=original_id)
    session.add(processed_group)

    session.commit()


def processed_stars_associations(original_stars: List[Star],
                                 processed_stars: List[Star]
                                 ) -> Iterator[associations.ProcessedStars]:
    for original_star, processed_star in zip(original_stars,
                                             processed_stars):
        yield associations.ProcessedStars(
                original_id=original_star.id,
                processed_id=processed_star.id)

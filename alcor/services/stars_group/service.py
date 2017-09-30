import uuid
from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import (Iterator,
                    List)

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star,
                          eliminations,
                          associations)
from alcor.models.star import set_radial_velocity_to_zero
from alcor.services.data_access import fetch_unprocessed_stars
from . import elimination


def process(*,
            group: Group,
            filtration_method: str,
            session: Session) -> None:
    stars = fetch_unprocessed_stars(group_id=group.id,
                                    session=session)
    stars_count = len(stars)

    eliminations_counter = Counter()
    if filtration_method in {'restricted', 'full'}:
        is_eliminated = partial(elimination.check,
                                eliminations_counter=eliminations_counter,
                                filtration_method=filtration_method)
        stars = list(filterfalse(is_eliminated, stars))

    counter = eliminations.StarsCounter(group_id=group.id,
                                        raw=stars_count,
                                        **eliminations_counter)
    session.add(counter)

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

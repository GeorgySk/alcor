import logging
from typing import List

from cassandra.cluster import Session

from alcor.models.group import Group
from alcor.models.star import Star
from alcor.services.data_access.reading import fetch


logger = logging.getLogger(__name__)


def delete_processed_data(session: Session) -> None:
    prepared_delete_star_statement = session.prepare(
        """
            DELETE FROM  alcor.stars WHERE id=?
        """
    )
    prepared_delete_group_statement = session.prepare(
        """
            DELETE FROM  alcor.groups WHERE id=?
        """
    )

    processed_groups = fetch_processed_groups(session=session)
    logger.debug(f'Going to delete information about {len(processed_groups)}'
                 f'processed groups')

    for group_count, group in enumerate(processed_groups):
        logger.debug(f'Group Nº {group_count}')
        stars = fetch_stars(group=group,
                            session=session)
        for star in stars:
            session.execute(prepared_delete_star_statement,
                            [star.id])
        logger.debug(f'Stars from group Nº {group_count} are deleted')
        logger.debug(f'Deleting group')
        session.execute(prepared_delete_group_statement,
                        [group.id])

    logger.debug(f'Finished deleting stars')


# TODO: I already have it, but I had circular dependencies so I left it here
def fetch_processed_groups(*,
                           session: Session
                           ) -> List[Group]:
    query = (Group.objects
             .filter(Group.processed == True)
             .limit(None))
    records = fetch(query=query,
                    session=session)
    return [Group(**record)
            for record in records]


# TODO: I already have it, but I had circular dependencies so I left it here
def fetch_stars(*,
                group: Group,
                session: Session
                ) -> List[Star]:
    logger.debug(f'Fetching stars')
    query = (Star.objects
             .filter(Star.group_id == group.id)
             .limit(None))
    logger.debug(f'Query for fetching stars is {query}')
    records = fetch(query=query,
                    session=session)
    return [Star(**record) for record in records]

from datetime import datetime
from typing import List

from cassandra.cluster import Session

from alcor.models import (Group,
                          Star)
from alcor.services.data_access import (fetch,
                                        model_update_statement)
from alcor.services.data_access.execution import execute_prepared_statement
from alcor.services.group_processing import process_stars_group


def run_processing(*,
                   filtration_method: str,
                   nullify_radial_velocity: bool,
                   luminosity_function: bool,
                   velocities_clouds: bool,
                   velocities_vs_magnitude: bool,
                   lepine_criterion: bool,
                   session: Session) -> None:
    groups = fetch_unprocessed_groups(session=session)
    for group in groups:
        stars = fetch_stars(group=group,
                            session=session)
        update_group(group,
                     session=session)
        process_stars_group(stars=stars,
                            group=group,
                            filtration_method=filtration_method,
                            nullify_radial_velocity=nullify_radial_velocity,
                            luminosity_function=luminosity_function,
                            velocities_clouds=velocities_clouds,
                            velocities_vs_magnitude=velocities_vs_magnitude,
                            lepine_criterion=lepine_criterion,
                            session=session)


def update_group(group: Group,
                 *,
                 session: Session) -> None:
    columns = [Group.processed, Group.updated_timestamp]
    where_clauses = [Group.id == group.id]
    statement = model_update_statement(Group,
                                       columns=columns,
                                       where_clauses=where_clauses,
                                       include_keyspace=False)
    execute_prepared_statement(statement=statement,
                               parameters_collection=[(True, datetime.now())],
                               session=session)


def fetch_unprocessed_groups(*,
                             session: Session
                             ) -> List[Group]:
    query = (Group.objects
             .filter(Group.processed == False)
             .limit(None))
    records = fetch(query=query,
                    session=session)
    return [Group(**record)
            for record in records]


def fetch_stars(*,
                group: Group,
                session: Session
                ) -> List[Star]:
    query = (Star.objects
             .filter(Star.group_id == group.id)
             .limit(None))
    records = fetch(query=query,
                    session=session)
    return [Star(**record) for record in records]

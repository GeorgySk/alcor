from datetime import datetime
from typing import (List,
                    Tuple)

from cassandra.cluster import Session as CassandraSession
from sqlalchemy.orm.session import Session

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.models import (Group,
                          Star)
from alcor.services.data_access import (fetch,
                                        model_update_statement)
from alcor.services.data_access.execution import execute_prepared_statement
from alcor.services.group_processing import process_stars_group


def run_processing(*,
                   filtration_method: str,
                   nullify_radial_velocity: bool,
                   w_luminosity_function: bool,
                   w_velocities_clouds: bool,
                   w_velocities_vs_magnitude: bool,
                   w_lepine_criterion: bool,
                   session: Session,
                   cassandra_session: CassandraSession) -> None:
    c_groups, groups = fetch_unprocessed_groups(
        cassandra_session=cassandra_session,
        session=session)
    for c_group, group in zip(c_groups, groups):
        c_stars, stars = fetch_stars(c_group=c_group,
                                     group=group,
                                     cassandra_session=cassandra_session,
                                     session=session)
        update_group(c_group,
                     group,
                     cassandra_session=cassandra_session,
                     session=session)
        process_stars_group(c_stars=c_stars,
                            stars=stars,
                            c_group=c_group,
                            group=group,
                            filtration_method=filtration_method,
                            nullify_radial_velocity=nullify_radial_velocity,
                            w_luminosity_function=w_luminosity_function,
                            w_velocities_clouds=w_velocities_clouds,
                            w_velocities_vs_magnitude=w_velocities_vs_magnitude,
                            w_lepine_criterion=w_lepine_criterion,
                            cassandra_session=cassandra_session,
                            session=session)


def update_group(c_group: CGroup,
                 group: Group,
                 *,
                 cassandra_session: CassandraSession,
                 session: Session) -> None:
    current_date_time = datetime.now()
    columns = [CGroup.processed, CGroup.updated_timestamp]
    where_clauses = [CGroup.id == c_group.id]
    statement = model_update_statement(CGroup,
                                       columns=columns,
                                       where_clauses=where_clauses,
                                       include_keyspace=False)
    execute_prepared_statement(statement=statement,
                               parameters_collection=[(True,
                                                       current_date_time)],
                               session=cassandra_session)

    query = (session.query(Group)
             .filter(Group.id == group.id))
    query.update({Group.processed: True,
                  Group.updated_timestamp: current_date_time})


def fetch_unprocessed_groups(*,
                             cassandra_session: CassandraSession,
                             session: Session
                             ) -> Tuple[List[CGroup],
                                        List[Group]]:
    query = (session.query(Group)
             .filter(Group.processed.is_(False)))
    c_query = (CGroup.objects
               .filter(CGroup.processed == False)
               .limit(None))
    records = fetch(query=c_query,
                    session=cassandra_session)
    return ([CGroup(**record)
             for record in records],
            query.all())


def fetch_stars(*,
                c_group: CGroup,
                group: Group,
                session: Session,
                cassandra_session: CassandraSession
                ) -> Tuple[List[CStar],
                           List[Star]]:
    query = (session.query(Star)
             .filter(Star.group_id == group.id))
    c_query = (CStar.objects
               .filter(CStar.group_id == c_group.id)
               .limit(None))
    records = fetch(query=c_query,
                    session=cassandra_session)
    return ([CStar(**record) for record in records],
            query.all())

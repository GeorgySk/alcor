from typing import (Set,
                    List)
from uuid import UUID

from cassandra.cluster import Session

from alcor.models import (Star,
                          simulation)
from alcor.services.data_access import fetch_query
from alcor.services.group_processing import process_stars_group


def run_processing(*,
                   filtration_method: str,
                   nullify_radial_velocity: bool,
                   luminosity_function: bool,
                   velocities_clouds: bool,
                   velocities_vs_magnitude: bool,
                   lepine_criterion: bool,
                   session: Session) -> None:
    groups_ids = fetch_groups_ids(session=session)
    for group_id in groups_ids:
        stars = fetch_stars(group_id=group_id,
                            session=session)
        process_stars_group(stars=stars,
                            group_id=group_id,
                            filtration_method=filtration_method,
                            nullify_radial_velocity=nullify_radial_velocity,
                            luminosity_function=luminosity_function,
                            velocities_clouds=velocities_clouds,
                            velocities_vs_magnitude=velocities_vs_magnitude,
                            lepine_criterion=lepine_criterion,
                            session=session)


def fetch_groups_ids(*,
                     session: Session
                     ) -> Set[UUID]:
    group_id_column = simulation.Parameter.group_id.column
    group_id_column_name = group_id_column.column_name
    query = (simulation.Parameter.objects
             .only([group_id_column_name])
             .limit(None))
    records = fetch_query(query=query,
                          session=session)
    return {record[group_id_column_name]
            for record in records}


def fetch_stars(*,
                group_id: UUID,
                session: Session
                ) -> List[Star]:
    query = (Star.objects
             .filter(Star.group_id == group_id)
             .limit(None))
    records = fetch_query(query=query,
                          session=session)
    return [Star(**record) for record in records]

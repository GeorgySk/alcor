import uuid
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models.group import Group
from alcor.models.luminosity_function import Point
from alcor.models.star import Star
from alcor.models.velocities.clouds import (Cloud,
                                            LepineCaseUVCloud,
                                            LepineCaseUWCloud,
                                            LepineCaseVWCloud)
from alcor.models.velocities_vs_magnitudes.bins import (Bin,
                                                        LepineCaseUBin,
                                                        LepineCaseVBin,
                                                        LepineCaseWBin)
from alcor.models.velocities_vs_magnitudes.clouds import (LepineCaseUCloud,
                                                          LepineCaseVCloud,
                                                          LepineCaseWCloud)


# TODO: implement only 1 function for all these
def fetch_stars(*,
                group: Group,
                session: Session
                ) -> List[Star]:
    query = (session.query(Star)
             .filter(Star.group_id == group.id))
    return query.all()


def fetch_unprocessed_groups(session: Session) -> List[Group]:
    query = (session.query(Group)
             .filter(Group.original_id.is_(None)))
    return query.all()


def fetch_last_groups(*,
                      count: int,
                      session: Session) -> List[Group]:
    query = (session.query(Group)
             .order_by(Group.updated_timestamp.desc())
             .limit(count))
    return query.all()


def fetch_group_by_id(*,
                      group_id: uuid.UUID,
                      session: Session) -> List[Group]:
    query = (session.query(Group)
             .filter(Group.id == group_id))
    return query.all()


def fetch_all_graph_points(session: Session) -> List[Point]:
    query = session.query(Point)
    return query.all()


def fetch_all_stars(session: Session) -> List[Star]:
    query = session.query(Star)
    return query.all()


def fetch_all_cloud_points(session: Session) -> List[Cloud]:
    query = session.query(Cloud)
    return query.all()


# TODO: this is the same as fetch_stars
def fetch_stars_by_group_id(*,
                            group_id: uuid.UUID,
                            session: Session) -> List[Star]:
    query = (session.query(Star)
             .filter(Star.group_id == group_id))
    return query.all()


def fetch_all_bins(session: Session) -> List[Bin]:
    query = session.query(Bin)
    return query.all()


def fetch_all_u_vs_mag_bins(session: Session) -> List[LepineCaseUBin]:
    query = session.query(LepineCaseUBin)
    return query.all()


def fetch_all_v_vs_mag_bins(session: Session) -> List[LepineCaseVBin]:
    query = session.query(LepineCaseVBin)
    return query.all()


def fetch_all_w_vs_mag_bins(session: Session) -> List[LepineCaseWBin]:
    query = session.query(LepineCaseWBin)
    return query.all()


def fetch_all_u_vs_mag_clouds(session: Session) -> List[LepineCaseUCloud]:
    query = session.query(LepineCaseUCloud)
    return query.all()


def fetch_all_v_vs_mag_clouds(session: Session) -> List[LepineCaseVCloud]:
    query = session.query(LepineCaseVCloud)
    return query.all()


def fetch_all_w_vs_mag_clouds(session: Session) -> List[LepineCaseWCloud]:
    query = session.query(LepineCaseWCloud)
    return query.all()


def fetch_all_lepine_case_uv_cloud_points(session: Session
                                          ) -> List[LepineCaseUVCloud]:
    query = session.query(LepineCaseUVCloud)
    return query.all()


def fetch_all_lepine_case_uw_cloud_points(session: Session
                                          ) -> List[LepineCaseUWCloud]:
    query = session.query(LepineCaseUWCloud)
    return query.all()


def fetch_all_lepine_case_vw_cloud_points(session: Session
                                          ) -> List[LepineCaseVWCloud]:
    query = session.query(LepineCaseVWCloud)
    return query.all()

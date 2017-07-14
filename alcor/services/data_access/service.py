from cassandra.cluster import Session

from alcor.services.data_access.deleting import delete_processed_data


def run_db_managing(with_deleting_processed_data: bool,
                    session: Session) -> None:
    if with_deleting_processed_data:
        delete_processed_data(session)

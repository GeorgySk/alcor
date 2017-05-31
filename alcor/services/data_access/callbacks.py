from threading import Event
from typing import Callable, Optional, List

from cassandra.cluster import ResponseFuture

from alcor.types import CallbackType, ResponseType, RecordType


class PagedResultHandler:
    """
    more info at
    http://datastax.github.io/python-driver/query_paging.html#handling-paged-results-with-callbacks
    """

    def __init__(self,
                 future: ResponseFuture,
                 callback: CallbackType,
                 event_type: Callable[[], Event] = Event):
        self.error = None
        self.callback = callback
        self.finished_event = event_type()
        self.future = future
        self.future.add_callbacks(
            callback=self.handle_page,
            errback=self.handle_error)

    def handle_page(self, records: ResponseType) -> None:
        self.callback(records)
        if self.future.has_more_pages:
            self.future.start_fetching_next_page()
        else:
            self.finished_event.set()

    def handle_error(self, exc: Exception) -> None:
        self.error = exc
        self.finished_event.set()


def empty_callback(records: ResponseType) -> None:
    pass


def add_callback(*,
                 future: ResponseFuture,
                 callback: CallbackType
                 ) -> Optional[List[RecordType]]:
    if callback is not None:
        handler = PagedResultHandler(future,
                                     callback=callback)
        handler.finished_event.wait()
        return

    return future.result()

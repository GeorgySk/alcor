from typing import Iterator


def file_reader(path: str) -> Iterator[str]:
    with open(path, mode='r') as file:
        yield from file

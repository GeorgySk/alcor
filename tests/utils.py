import inspect
import os
from typing import (Any,
                    Callable,
                    Iterable,
                    Hashable,
                    Dict,
                    List)

from hypothesis import (Verbosity,
                        find,
                        settings)
from hypothesis.searchstrategy import SearchStrategy


def example(strategy: SearchStrategy) -> Any:
    return find(specifier=strategy,
                condition=lambda x: True,
                settings=settings(max_shrinks=0,
                                  max_iterations=10000,
                                  database=None,
                                  verbosity=Verbosity.quiet))


def double_star_map(function: Callable[..., Any],
                    kwargs: Dict[str, Any]) -> Any:
    return function(**kwargs)


def initializer_parameters(cls: type) -> List[str]:
    initializer_signature = inspect.signature(cls.__init__)
    parameters = dict(initializer_signature.parameters)
    parameters.pop('self')
    return list(parameters.keys())


def sub_dict(dictionary: Dict[Hashable, Any],
             *,
             keys: Iterable[Hashable]) -> Dict[Hashable, Any]:
    return {key: dictionary[key]
            for key in keys}


def fort_files_by_metallicities_lengths(fort_links: Dict[int, Iterable[int]],
                                        *,
                                        base_dir: str) -> List[int]:
    for fort_links_range in fort_links.values():
        yield from fort_files_lengths(fort_links=fort_links_range,
                                      base_dir=base_dir)


def fort_files_lengths(fort_links: Iterable[int],
                       *,
                       base_dir: str) -> List[int]:
    fort_link_dirs = [os.path.join(base_dir,
                                   'fort_files/fort.' + str(fort_link))
                      for fort_link in fort_links]
    yield from map(file_lines_count, fort_link_dirs)


def file_lines_count(file_path: str) -> int:
    with open(file_path) as file:
        return sum(1 for _ in file)

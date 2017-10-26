import inspect
import os
from typing import (Any,
                    Callable,
                    Iterable,
                    Iterator,
                    Hashable,
                    Dict,
                    List,)

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


def fort_files_by_metallicities_lengths(fort_links: Dict[int, Iterator],
                                        *,
                                        base_dir: str) -> List[int]:
    lengths = []
    for metallicity, fort_links_range in fort_links.items():
        lengths.extend(fort_files_lengths(fort_links=fort_links_range,
                                          base_dir=base_dir))
    return lengths


def fort_files_lengths(fort_links: Iterator,
                       *,
                       base_dir: str) -> List[int]:
    lengths = []
    fort_link_dirs = [os.path.join(base_dir,
                                   'fort_files/fort.' + str(fort_link))
                      for fort_link in fort_links]
    for fort_link_dir in fort_link_dirs:
        lengths.append(sum(1 for line in open(fort_link_dir)))

    return lengths

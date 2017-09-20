import sys
from decimal import Decimal
from functools import partial
from typing import (Iterator,
                    Tuple)

from hypothesis import strategies
from hypothesis.searchstrategy import SearchStrategy
from hypothesis_sqlalchemy import tables
from sqlalchemy import Float

from alcor.models import Star


def stars_factory(nullable: bool) -> SearchStrategy:
    table = Star.__table__
    fixed_columns_values = dict(fixed_stars_columns_values(nullable))
    records = tables.records.factory(table, **fixed_columns_values)

    columns_names = [column.name
                     for column in table.columns]
    return (records
            .map(partial(zip, columns_names))
            .map(dict)
            .map(Star.deserialize))


def fixed_stars_columns_values(nullable: bool
                               ) -> Iterator[Tuple[str, SearchStrategy]]:
    table = Star.__table__
    positive_decimals = strategies.decimals(
            min_value=Decimal(sys.float_info.epsilon),
            allow_nan=False,
            allow_infinity=False)
    if nullable:
        positive_decimals |= strategies.none()
    else:
        decimals = strategies.decimals(allow_nan=False,
                                       allow_infinity=False)
        for column in table.columns:
            if isinstance(column.type, Float):
                yield column.name, decimals
    yield Star.proper_motion.name, positive_decimals


defined_stars = stars_factory(nullable=False)
undefined_stars = stars_factory(nullable=True)

defined_stars_lists = strategies.lists(defined_stars)
undefined_stars_lists = strategies.lists(undefined_stars)

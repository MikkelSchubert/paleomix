# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import bisect
import random
from collections.abc import Iterable, Iterator
from typing import TypeVar

T = TypeVar("T")


def weighted_sampling(
    choices: Iterable[T],
    weights: Iterable[float],
    rng: random.Random | None = None,
) -> Iterator[T]:
    if rng is None:
        rng = random.Random()

    choices = list(choices)
    weights = list(weights)
    if not weights or (len(weights) != len(choices)):
        raise ValueError(
            "Choices and probabilities must be non-empty lists "
            f"of identical length, not lengths {len(choices)} and {len(weights)}"
        )

    total: float = 0
    totals: list[float] = []
    for index, weight in enumerate(weights, start=1):
        if weight <= 0:
            raise ValueError(
                f"Probabilities must be > 0, not {weight!r} for weight {index}"
            )
        total += weight
        totals.append(total)

    while True:
        rand = rng.random() * total
        index = bisect.bisect_right(totals, rand)
        yield choices[index]


def reservoir_sampling(
    items: Iterable[T],
    downsample_to: int,
    rng: random.Random | None = None,
) -> list[T]:
    if rng is None:
        rng = random.Random()

    if not isinstance(downsample_to, int):
        raise TypeError(f"downsample_to must be an int, not {downsample_to!r}")
    elif downsample_to < 0:
        raise ValueError(f"downsample_to must be >= 0, not {downsample_to}")

    reservoir: list[T] = []
    for index, item in enumerate(items):
        if index >= downsample_to:
            index = rng.randint(0, index)
            if index < downsample_to:
                reservoir[index] = item
        else:
            reservoir.append(item)
    return reservoir

#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
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
            "of identical length, not lengths %i and %i" % (len(choices), len(weights))
        )

    total: float = 0
    totals: list[float] = []
    for index, weight in enumerate(weights, start=1):
        if weight <= 0:
            raise ValueError(
                "Probabilities must be > 0, not %r for weight %i" % (weight, index)
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

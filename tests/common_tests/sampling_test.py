# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

from collections.abc import Iterable
from unittest.mock import Mock

import pytest

from paleomix.common import sampling

###############################################################################
###############################################################################
# weighted_sampling

_SELECT_BY_VALUE = (
    (0.00000, "a"),
    (0.16666, "a"),  # < 1/6
    (1 / 6.0, "b"),
    (0.49999, "b"),  # < 3/6
    (3 / 6.0, "c"),
    (0.99999, "c"),
)


@pytest.mark.parametrize(("value", "expectation"), _SELECT_BY_VALUE)
def test_weighted_sampling__select_by_weight(value: float, expectation: str) -> None:
    choices = "abc"
    weights = (1, 2, 3)
    rng = Mock(random=lambda: value)
    iterator = sampling.weighted_sampling(choices, weights, rng)
    assert next(iterator) == expectation


@pytest.mark.parametrize(("choices", "weights"), [([], []), ([], [1, 2]), ([1, 2], [])])
def test_weighted_sampling__empty_lists_raises(
    choices: list[int],
    weights: list[float],
) -> None:
    iterator = sampling.weighted_sampling(choices, weights)
    with pytest.raises(ValueError, match="Choices and probabilities must be non-empty"):
        iterator.__next__()


_MISMATCHED_LENGTH_INPUTS = (
    ([0, 1], [1, 2, 3]),
    ([0, 1, 2], [1, 2]),
    (iter([0, 1]), [1, 2, 3]),
    ([0, 1], iter([1, 2, 3])),
    (iter([0, 1]), iter([1, 2, 3])),
)


@pytest.mark.parametrize(("choices", "weights"), _MISMATCHED_LENGTH_INPUTS)
def test_weighted_sampling__different_length_input_raises_value_error(
    choices: Iterable[int],
    weights: list[float],
) -> None:
    iterator = sampling.weighted_sampling(choices, weights)
    with pytest.raises(ValueError, match="Choices and probabilities must be non-empty"):
        iterator.__next__()


def test_weighted_sampling__negative_weight_value_error() -> None:
    choices = [0, 1, 2]
    weights = [1, -2, 3]
    iterator = sampling.weighted_sampling(choices, weights)
    with pytest.raises(ValueError, match="Probabilities must be > 0"):
        iterator.__next__()


def test_weighted_sampling__zero_weight_raises_value_error() -> None:
    choices = [0, 1, 2]
    weights = [1, 0, 3]
    iterator = sampling.weighted_sampling(choices, weights)
    with pytest.raises(ValueError, match="Probabilities must be > 0"):
        iterator.__next__()


def test_weighted_sampling__non_numerical_weight_raises_type_error() -> None:
    choices = [0, 1, 2]
    weights = [1, "foo", 3]
    iterator = sampling.weighted_sampling(
        choices,
        weights,  # pyright: ignore[reportGeneralTypeIssues]
    )
    with pytest.raises(TypeError):
        iterator.__next__()


###############################################################################
###############################################################################
# reservoir_sampling


def test_reservoir_sampling__select_first_item() -> None:
    def randint(_min: int, _max: int) -> int:
        return 1

    rng = Mock(randint=randint)
    values = [1, 2]
    result = sampling.reservoir_sampling(values, 1, rng)
    assert result == [1]


def test_reservoir_sampling__select_second_item() -> None:
    def randint(_min: int, _max: int) -> int:
        return 0

    rng = Mock(randint=randint)
    values = [1, 2]
    result = sampling.reservoir_sampling(values, 1, rng)
    assert result == [2]


def test_reservoir_sampling__upsample_equals_input() -> None:
    result = sampling.reservoir_sampling(list(range(5)), 10)
    assert result == list(range(5))


def test_reservoir_sampling__downsample_to_zero() -> None:
    result = sampling.reservoir_sampling(list(range(5)), 0)
    assert result == []


def test_reservoir_sampling__downsample_to_negative_raises_value_error() -> None:
    with pytest.raises(ValueError, match="downsample_to must be >= 0"):
        sampling.reservoir_sampling(list(range(5)), -1)


def test_reservoir_sampling__downsample_to_float_raises_type_error() -> None:
    with pytest.raises(TypeError):
        sampling.reservoir_sampling(
            list(range(5)),
            1.0,  # pyright: ignore[reportArgumentType]
        )


def test_reservoir_sampling__downsample_to_non_number_raises_type_error() -> None:
    with pytest.raises(TypeError):
        sampling.reservoir_sampling(
            list(range(5)),
            "Eh?",  # pyright: ignore[reportArgumentType]
        )

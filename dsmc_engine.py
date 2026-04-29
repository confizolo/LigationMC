"""Event-driven DSMC gelation engine using Gillespie SSA."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np

from polymer_utils import cyclisation_rate, smoluchowski_kernel


@dataclass
class MergeEvent:
    time: float
    length_i: int
    length_j: int
    new_length: int


@dataclass
class CyclisationEvent:
    time: float
    linear_length: int
    ring_length: int
    links_formed: int = 0
    linked_ring_ids: list[int] | None = None
    ring_id: int | None = None


Event = Union[MergeEvent, CyclisationEvent]


class DSMCEngine:
    """Stochastic DSMC engine for linear-polymer gelation."""

    def __init__(
        self,
        linear_lengths: list[int],
        k1: float,
        k2: float,
        alpha: float = 1.0,
        nu: float = 0.5,
        seed: int | None = None,
        rng: np.random.Generator | None = None,
    ) -> None:
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.alpha = float(alpha)
        self.nu = float(nu)
        self.time = 0.0

        self._linear_population: dict[int, int] = {}
        for length in linear_lengths:
            self._linear_population[length] = self._linear_population.get(length, 0) + 1

        self._ring_population: list[int] = []
        self._rng = rng if rng is not None else np.random.default_rng(seed)

    @property
    def linear_population(self) -> dict[int, int]:
        return dict(self._linear_population)

    @property
    def ring_population(self) -> list[int]:
        return list(self._ring_population)

    @property
    def n_linear(self) -> int:
        return sum(self._linear_population.values())

    def _reaction_channels(self) -> tuple[list[tuple[str, int, int]], np.ndarray]:
        lengths = sorted(self._linear_population)
        channels: list[tuple[str, int, int]] = []
        propensities: list[float] = []

        for idx, i in enumerate(lengths):
            ni = self._linear_population[i]
            if ni <= 0:
                continue

            a_cyc = ni * cyclisation_rate(i, self.k2, self.nu)
            if a_cyc > 0.0:
                channels.append(("cyc", i, i))
                propensities.append(a_cyc)

            for j in lengths[idx:]:
                nj = self._linear_population[j]
                if nj <= 0:
                    continue
                if i == j and ni < 2:
                    continue

                factor = ni * (nj - (1 if i == j else 0))
                if factor <= 0:
                    continue
                a_merge = self.k1 * factor * smoluchowski_kernel(i, j, self.alpha, self.nu)
                if a_merge > 0.0:
                    channels.append(("merge", i, j))
                    propensities.append(a_merge)

        return channels, np.asarray(propensities, dtype=float)

    def step(self) -> Event:
        if self.n_linear <= 0:
            raise RuntimeError("No linear polymers left to evolve.")

        channels, prop = self._reaction_channels()
        if prop.size == 0:
            raise RuntimeError("No available reaction channels while linears remain.")

        a0 = float(prop.sum())
        tau = float(self._rng.exponential(1.0 / a0))
        self.time += tau

        choice = float(self._rng.random() * a0)
        cumulative = 0.0
        picked = 0
        for idx, a in enumerate(prop):
            cumulative += float(a)
            if choice <= cumulative:
                picked = idx
                break

        channel, i, j = channels[picked]

        if channel == "merge":
            self._linear_population[i] -= 1
            self._linear_population[j] -= 1
            if self._linear_population[i] == 0:
                del self._linear_population[i]
            if j in self._linear_population and self._linear_population[j] == 0:
                del self._linear_population[j]

            new_length = i + j
            self._linear_population[new_length] = self._linear_population.get(new_length, 0) + 1
            return MergeEvent(time=self.time, length_i=i, length_j=j, new_length=new_length)

        self._linear_population[i] -= 1
        if self._linear_population[i] == 0:
            del self._linear_population[i]

        self._ring_population.append(i)
        return CyclisationEvent(time=self.time, linear_length=i, ring_length=i, linked_ring_ids=[])

    def run_until_exhausted(self, max_steps: int = 50000) -> list[Event]:
        events: list[Event] = []
        steps = 0
        while self.n_linear > 0:
            if steps >= max_steps:
                raise RuntimeError(f"Exceeded max_steps={max_steps} before exhausting linears.")
            events.append(self.step())
            steps += 1
        return events

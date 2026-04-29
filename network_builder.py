"""Graph-based ring-network construction for LigMC."""

from __future__ import annotations

from collections import Counter
import math
from dataclasses import replace

import networkx as nx
import numpy as np

from dsmc_engine import CyclisationEvent
from polymer_utils import valence_model


class NetworkBuilder:
    """Build and query the evolving graph of topologically linked rings."""

    def __init__(
        self,
        initial_ring_lengths: list[int],
        seed: int | None = None,
        rng: np.random.Generator | None = None,
    ) -> None:
        self._graph = nx.Graph()
        self._ring_lengths: dict[int, int] = {}
        self._rng = rng if rng is not None else np.random.default_rng(seed)

        for length in initial_ring_lengths:
            self._add_ring(length)

    def _add_ring(self, length: int) -> int:
        ring_id = self._graph.number_of_nodes()
        self._graph.add_node(ring_id, length=int(length))
        self._ring_lengths[ring_id] = int(length)
        return ring_id

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @property
    def ring_lengths(self) -> dict[int, int]:
        return dict(self._ring_lengths)

    def process_cyclisation(
        self,
        event: CyclisationEvent,
        A: float,
        box_volume: float = 1.0,
    ) -> CyclisationEvent:
        new_ring_id = self._add_ring(event.ring_length)
        l_cyc = int(event.ring_length)

        pool = [node for node in self._graph.nodes if node != new_ring_id]
        n_total = len(pool)
        if not pool:
            return replace(event, links_formed=0, linked_ring_ids=[], ring_id=new_ring_id)

        # Per-target linking probability derived from Poisson(mu) with
        # mu = A * nring_target * l_cyc / box_volume. Probability of at least
        # one event is p = 1 - exp(-mu). Draw independently for each existing
        # ring and add a single edge when the draw is successful.
        targets: list[int] = []
        for target in pool:
            nring_target = int(self._ring_lengths.get(int(target), 0))
            if nring_target <= 0:
                continue
            mu = float(A) * float(nring_target) * float(l_cyc) / float(box_volume)
            if mu <= 0.0:
                continue
            p = 1.0 - math.exp(-mu)
            if self._rng.random() < p:
                self._graph.add_edge(new_ring_id, int(target))
                targets.append(int(target))

        return replace(
            event,
            links_formed=len(targets),
            linked_ring_ids=[int(t) for t in targets],
            ring_id=new_ring_id,
        )

    def largest_component_fraction(self) -> float:
        if self._graph.number_of_nodes() == 0:
            return 0.0
        largest = max(nx.connected_components(self._graph), key=len)
        return len(largest) / self._graph.number_of_nodes()

    def degree_distribution(self) -> dict[int, int]:
        degrees = [self._graph.degree(node) for node in self._graph.nodes]
        return dict(Counter(degrees))

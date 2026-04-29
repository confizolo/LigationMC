"""Shared polymer-physics helpers for LigMC simulations."""

from __future__ import annotations

import math
import os

RESULTS_DIR = "/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc"
os.makedirs(RESULTS_DIR, exist_ok=True)


def calculate_valence(nlin: int, nring: int, blin: int = 64, bring: int = 512) -> float:
    """Return the mean of the Poisson valence distribution used for linking."""
    if nlin <= 0:
        return 0.0
    return float((nlin / blin) * (nring / bring) ** -0.6)


def calculate_polymer_numbers(
    ccsr: float,
    ccsl: float,
    L: float = 200,
    Rgr_base: float = 9.9,
    Rgl_base: float = 14.0,
    Nr: int = 1024,
    Nl: int = 128,
) -> tuple[int, int]:
    """Compute (Mr, Ml) using overlap-concentration scaling from legacy code."""
    rgr_scaled = Rgr_base * math.sqrt(Nr / 128.0)
    rgl_scaled = Rgl_base * math.sqrt(Nl / 128.0)

    v_box = L**3
    v_polymer_ring = (4.0 / 3.0) * math.pi * (rgr_scaled**3)
    v_polymer_linear = (4.0 / 3.0) * math.pi * (rgl_scaled**3)

    mrc = v_box / v_polymer_ring
    mlc = v_box / v_polymer_linear

    mr = math.floor(mrc * ccsr)
    ml = math.floor(mlc * ccsl)
    return mr, ml


def smoluchowski_kernel(i: int, j: int, alpha: float = 1.0, nu: float = 0.5) -> float:
    """Generalized Smoluchowski kernel from De Gennes scaling."""
    return (i ** (-alpha) + j ** (-alpha)) * (i**nu + j**nu)


def cyclisation_rate(length: int, k2: float, nu: float = 0.5) -> float:
    """Return the single-polymer cyclisation propensity contribution."""
    return k2 * (length ** (-4.0 * nu))


def valence_model(l_cyc: int, n_total: int, A: float, box_volume: float = 1.0) -> float:
    """Concentration-driven Poisson mean for topological linking.

    B is hardcoded to 1.0, so:
    lambda = A * ((N_total * l_cyc) / V_box)
    """
    if l_cyc <= 0 or n_total <= 0:
        return 0.0
    if box_volume <= 0.0:
        return 0.0
    concentration_proxy = (n_total * l_cyc) / box_volume
    return float(A * concentration_proxy)

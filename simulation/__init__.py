"""LigMC simulation core — DSMC engine, network builder, and fitting tools."""

from .dsmc_engine import DSMCEngine, MergeEvent, CyclisationEvent, Event
from .network_builder import NetworkBuilder
from .polymer_utils import (
    RESULTS_DIR,
    calculate_polymer_numbers,
    smoluchowski_kernel,
    cyclisation_rate,
    valence_model,
)
from .analysis import (
    save_trial_results,
    save_results_all,
    save_event_timeline,
    save_json,
    plot_gelation_curves,
    plot_degree_distributions,
    plot_stages_vs_nlin_nring,
)

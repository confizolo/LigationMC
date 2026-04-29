"""LigMC simulation core — DSMC engine, network builder, and fitting tools."""

from simulation.dsmc_engine import DSMCEngine, MergeEvent, CyclisationEvent, Event
from simulation.network_builder import NetworkBuilder
from simulation.polymer_utils import (
    RESULTS_DIR,
    calculate_polymer_numbers,
    smoluchowski_kernel,
    cyclisation_rate,
    valence_model,
)
from simulation.analysis import (
    save_trial_results,
    save_results_all,
    save_event_timeline,
    save_json,
    plot_gelation_curves,
    plot_degree_distributions,
    plot_stages_vs_nlin_nring,
)

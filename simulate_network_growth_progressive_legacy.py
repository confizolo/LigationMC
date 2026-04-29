
import os
RESULTS_DIR = '/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels/results/'
os.makedirs(RESULTS_DIR, exist_ok=True)
#!/usr/bin/env python3
"""
Simulate network growth starting from N0 isolated nodes that cannot be linked among themselves.
At each timestep a new node is added that links to f existing nodes (uniform random, without replacement).
Measure the number of added nodes required until the largest connected component reaches >=50% of the total network size.

Usage:
    python simulate_network_growth.py

Parameters can be adjusted at the bottom of the file or passed via command-line arguments.
"""

import random
import argparse
import numpy as np
import networkx as nx
from tqdm import trange
from mpi4py import MPI
import statistics
import sys
import matplotlib.pyplot as plt
import sys, os
import math
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def calculate_valence(nlin, nring, blin=64, bring=512):
    """Calculate the valence of a node based on its type and neighborhood."""
    if nlin > 0:
        return (nlin / blin) * (nring / bring) ** -0.6
    return 0

def calculate_polymer_numbers(
    ccsr: float,
    ccsl: float,
    L: float = 200,
    Rgr_base: float = 9.9,
    Rgl_base: float = 14.0,
    Nr: int = 1024,
    Nl: int = 128
) -> tuple[int, int]:
    """
    Calculates the number of ring and linear polymers based on notebook logic.
    This function replicates the calculations from the provided Mathematica notebook
    to determine the number of polymers (Mr, Ml) for a given simulation box
    and concentration scalars.
    Args:
        ccsr: The "chosen concentration scalar" for ring polymers.
               (e.g., 5, meaning 5 times the overlap concentration) 
        ccsl: The "chosen concentration scalar" for linear polymers.
               (e.g., 0.05, meaning 0.05 times the overlap concentration) 
        L: The box size.
        Rgr_base: The radius of gyration for a ring polymer of length 128.
        Rgl_base: The radius of gyration for a linear polymer of length 128.
        Nr: The polymer length (monomer count) for rings.
        Nl: The polymer length (monomer count) for linear polymers.
    Returns:
        A tuple containing (Mr, Ml):
        Mr (int): The final number of ring polymers (floored).
        Ml (int): The final number of linear polymers (floored).
    """
    # --- 1. Calculate the scaled radius of gyration for both types ---
    # This scales the base radius (for length 128) to the target length (Nr or Nl)
    # The notebook uses Sqrt[N/128], assuming Rg ~ N^0.5
    Rgr_scaled = Rgr_base * math.sqrt(Nr / 128.0)  # [cite: 7]
    Rgl_scaled = Rgl_base * math.sqrt(Nl / 128.0)  # [cite: 8]
    # --- 2. Calculate the volume of the box ---
    V_box = L**3  # [cite: 7, 8]
    # --- 3. Calculate the effective volume of a single polymer coil ---
    V_polymer_ring = (4.0 / 3.0) * math.pi * (Rgr_scaled**3)  # [cite: 7]
    V_polymer_linear = (4.0 / 3.0) * math.pi * (Rgl_scaled**3)  # [cite: 8]
    # --- 4. Calculate the number of polymers at overlap concentration (C*) ---
    # This is the number of polymers that would fill the box at C*
    Mrc = V_box / V_polymer_ring  # [cite: 7, 10]
    Mlc = V_box / V_polymer_linear  # [cite: 8, 11]
    # --- 5. Calculate the final number of polymers based on chosen concentration ---
    # This uses the user-provided scalars (ccsr, ccsl)
    Mr = math.floor(Mrc * ccsr)  # [cite: 12, 14]
    Ml = math.floor(Mlc * ccsl)  # [cite: 12]
    return Mr, Ml

def run_trial(mring, nring, mlin, nlin, max_steps=100000, seed=None):
    """Run one growth trial.

    N0: initial isolated nodes (labels 0..N0-1). They remain unlinked among themselves.
    f: number of links the new node forms to existing nodes when added.
    max_steps: maximum number of added nodes to attempt.
    Returns the number of added nodes needed to reach >=50% largest component, or None if not reached.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    G = nx.Graph()
    # Seed the graph with the isolated ring nodes (never connect among themselves)
    G.add_nodes_from(np.arange(mring))

    total_nodes = mring
    target_fraction = 0.5

    # Track sizes to allow future mass-based metrics if needed
    # Store monomer counts so that mass-based criteria can be revisited later
    node_sizes = np.full(total_nodes, nring)
    
    largest_comp_size = 1 if mring > 0 else 0

    if largest_comp_size >= target_fraction * total_nodes:
        return 0,G

    for step in range(1, max_steps + 1):
        for i in range(mlin):
            new_node = total_nodes
            
            G.add_node(new_node)
            node_sizes = np.append(node_sizes, nlin)

            # Only ring nodes are eligible attachment targets in this toy model
            existing_nodes = list(np.arange(mring).astype(int).tolist())

            # Draw a mean valence-based target count while respecting availability
            k = min(calculate_valence(nlin, nring), len(existing_nodes))
            num_targets = round(min(np.random.poisson(k), len(existing_nodes)))
            
            # If num_targets > 0, sample without replacement
            if num_targets > 0:
                targets = random.sample(existing_nodes, num_targets)
            else:
                targets = []

            for t in targets:
                G.add_edge(new_node, t)
    
            total_nodes += 1

        # # compute size of largest connected component
        # largest_comp_size = len(max(nx.connected_components(G), key=len))

        # # retrieve indices of nodes in the largest connected component
        # if G.number_of_nodes() == 0:
        #     largest_comp_indices = np.array([], dtype=int)
        # else:
        #     largest_cc = max(nx.connected_components(G), key=len)
        #     largest_comp_indices = np.array(sorted(int(n) for n in largest_cc), dtype=int)
        
        # largest_comp_size = node_sizes[largest_comp_indices].sum() if len(largest_comp_indices) > 0 else 0

        # if largest_comp_size >= target_fraction * node_sizes.sum():
        #     # Return number of added nodes to reach target (i.e., step)
        #     return step, G
            
        # check whether all original nodes (0..mring-1) are in the largest connected component
        if G.number_of_nodes() == 0:
            largest_cc = set()
        else:
            largest_cc = max(nx.connected_components(G), key=len)
            # print(largest_cc)

        if mring == 0:
            return step, G

        if set(range(mring)).issubset(largest_cc):
            return step, G
        
    return None, G

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate network growth until largest component >=50%')
    parser.add_argument('--mring', type=int, default=222, help='Initial number of isolated nodes')
    parser.add_argument('--nring', type=int, default=256, help='Initial number of isolated nodes')
    parser.add_argument('--mlin', type=int, default=2, help='Initial number of isolated nodes')
    parser.add_argument('--nlin', type=int, default=128, help='Initial number of isolated nodes')

    parser.add_argument('--trials', type=int, default=200, help='Number of independent trials')
    parser.add_argument('--max_steps', type=int, default=50000, help='Maximum number of added nodes per trial')

    args = parser.parse_args()

    results = {}

    # Sweep a grid of chain lengths to mimic the simulation matrix from TAPLab
    nrings = np.arange(256, 1040, 16)
    nlins = np.arange(64, 224, 16)

    systems = []
    for nring in nrings:  
        for nlin in nlins:  
            mring, mlin = calculate_polymer_numbers(5, 0.05, Nr=nring, Nl=nlin)
            systems.append(f"80 {mring} {nring} {mlin} {nlin}")

    results_dir = '/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels/python_simpl_model/sim_network_growth'
    os.makedirs(results_dir, exist_ok=True)

    for system in systems:
        if rank == 0:
            print(f"Processing system: {system}")
        sets = np.array(system.split(" ")).astype(int)

        L = sets[0]
        mring = sets[1]
        nring = sets[2]
        mlin = sets[3]
        nlin = sets[4]

        # Cache per-system trials to enable incremental runs/restarts
        fname = os.path.join(results_dir, f"results_L{L}_mring{mring}_nring{nring}_mlin{mlin}_nlin{nlin}.pkl")

        # load existing results if present
        if os.path.exists(fname):
            with open(fname, 'rb') as f:
                existing = pickle.load(f)
            results[system] = existing
            results_local = []
            start_from = len(existing)
        else:
            results[system] = []
            results_local = []
            start_from = 0

        # Only run the missing trials so that long sweeps can resume quickly
        if start_from < args.trials:
            needed = args.trials - start_from

            perrank = needed // size
            remainder = needed % size

            results_global = np.zeros(needed, dtype=np.uint64)

            comm.Barrier()

            begin = rank * perrank + min(rank, remainder)
            end = begin + perrank + (1 if rank < remainder else 0)
            for i in range(begin, end):
                res, G = run_trial(mring, nring, mlin, nlin, max_steps=args.max_steps, seed=None)
                results_local.append(res)
                # incremental save after each trial to avoid data loss on interruption
            
            results_local = np.array(results_local, dtype=np.uint64)

            comm.Barrier()

            counts = np.array(comm.allgather(len(results_local)))
            displs = np.insert(np.cumsum(counts), 0, 0)[0:-1]

            comm.Allgatherv(results_local, [results_global, counts, displs, MPI.UNSIGNED_LONG])

            results[system].extend(results_global.tolist())
            
            if rank == 0:  
                with open(fname, 'wb') as f:
                    pickle.dump(results[system], f)
        else:
            # already have enough trials saved
            pass

    if rank == 0:
        # Save a consolidated dictionary for quick downstream analysis
        with open(os.path.join(results_dir, 'results_all.pkl'), 'wb') as f:
            pickle.dump(results, f)


        nlin_list = []
        mlin_list = []
        nring_list = []
        mring_list = []
        means = []
        fail_fractions = []

        # Collapse trial results into summary statistics for plotting
        for system, reslist in results.items():
            parts = system.split()
            # system format: "L mring nring mlin nlin"
            mring = int(parts[1])
            nring = int(parts[2])
            mlin = int(parts[3])
            nlin = int(parts[4])

            vals = np.array([r if r is not None else np.nan for r in reslist], dtype=float)
            med = np.nanmean(vals) if np.any(~np.isnan(vals)) else np.nan
            fail_frac = np.mean(np.isnan(vals))

            nlin_list.append(nlin)
            nring_list.append(nring)
            mlin_list.append(mlin)
            mring_list.append(mring)
            means.append(med)
            fail_fractions.append(fail_frac)

        nlin_arr = np.array(nlin_list)
        mlin_arr = np.array(mlin_list)
        nring_arr = np.array(nring_list)
        mring_arr = np.array(mring_list)
        mean_arr = np.array(means)
        valence_arr = np.array([calculate_valence(nlin, nring) for nlin, nring in zip(nlin_arr, nring_arr)])

        # Basic scatter: color = median stages (NaNs will appear blank)
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))

        sc = axes[0,0].scatter(valence_arr*mlin_arr, mean_arr, s=70, edgecolor='k')

        slopes = []

        # Fit per-nlin power laws to see if valence*count controls the stages
        for nlin in nlins:
            mask_nlin = (nlin_arr == nlin) 

            x_fit = np.log(valence_arr[mask_nlin]*mlin_arr[mask_nlin])
            y_fit = np.log(mean_arr[mask_nlin])
            coeffs = np.polyfit(x_fit, y_fit, 1)
            slopes.append((nlin, coeffs[0]))

            fit_y = np.exp(np.polyval(coeffs, x_fit))
            axes[0,0].plot(valence_arr[mask_nlin]*mlin_arr[mask_nlin], fit_y, '--', alpha=0.7, label=f'nlin={nlin} (slope={coeffs[0]:.2f})')


        axes[0,0].legend()
        axes[0,0].set_xscale('log')
        axes[0,0].set_yscale('log')

        axes[0,1].plot([nlin for nlin, slope in slopes], [slope for nlin, slope in slopes], 'o-')
        axes[0,1].set_xlabel('nlin')
        axes[0,1].set_ylabel('slope of log(stages)+1 vs log(nlin)')
        axes[0,1].set_yscale('log')
        axes[0,1].set_xscale('log')

        fit_coeffs = np.polyfit(np.log([nlin for nlin, slope in slopes]), np.log([slope for nlin, slope in slopes]), 1)
        fit_slope = np.exp(np.polyval(fit_coeffs, np.log([nlin for nlin, slope in slopes])))
        axes[0,1].plot([nlin for nlin, slope in slopes], fit_slope, 'r--', label=f'Power law fit (exp={fit_coeffs[0]:.2f})')
        axes[0,1].legend()

        # color by nlin_arr and add colorbar
        sc.set_array(nlin_arr)
        cb0 = fig.colorbar(sc, ax=axes[0,0])
        cb0.set_label('nlin')
        axes[0,0].set_xlabel('nlin*(nring^(-0.6)')
        axes[0,0].set_ylabel('Stages')
        # filter out NaN medians
        valid = ~np.isnan(mean_arr)

        # Median vs nlin (sorted) colored by nring
        order_nlin = np.argsort(nlin_arr[valid])
        x1 = nlin_arr[valid][order_nlin]
        y1 = mean_arr[valid][order_nlin]
        c1 = nring_arr[valid][order_nlin]
        sc1 = axes[1,0].scatter(x1, y1, c=c1, cmap='plasma', s=70, edgecolor='k')
        axes[1,0].set_xlabel('nlin')
        axes[1,0].set_ylabel('mean stages')
        axes[1,0].set_title('Mean stages vs nlin (color = nring)')
        cb1 = fig.colorbar(sc1, ax=axes[1,0])
        cb1.set_label('nring')

        # Median vs nring (sorted) colored by nlin
        order_nring = np.argsort(nring_arr[valid])
        mask = nlin_arr[valid][order_nring] > 0
        x2 = nring_arr[valid][order_nring][mask]
        y2 = mean_arr[valid][order_nring][mask]
        c2 = nlin_arr[valid][order_nring][mask]
        sc2 = axes[1,1].scatter(x2, y2, c=c2, cmap='viridis', s=70, edgecolor='k')
        axes[1,1].set_xlabel('nring')
        axes[1,1].set_ylabel('mean stages')
        axes[1,1].set_title('Mean stages vs nring (color = nlin)')
        axes[1,1].set_xscale('log')
        axes[1,1].set_yscale('log')

        for nlin in nlins:
            mask_nlin = (nlin_arr == nlin) & valid
            if np.sum(mask_nlin) < 2:
                continue
            x_fit = np.log(nring_arr[mask_nlin])
            y_fit = np.log(mean_arr[mask_nlin])
            coeffs = np.polyfit(x_fit, y_fit, 1)
            fit_y = np.exp(np.polyval(coeffs, x_fit))
            axes[1,1].plot(nring_arr[mask_nlin], fit_y, '--', alpha=0.7, label=f'nlin={nlin} (slope={coeffs[0]:.2f})')
        axes[1,1].legend()
        cb2 = fig.colorbar(sc2, ax=axes[1,1])
        cb2.set_label('nlin')

        print(fit_coeffs)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'stages_vs_nlin_nring.png'), dpi=150)
        plt.show()

        # Plot degree distribution for each system
        fig, axes = plt.subplots(4, 5, figsize=(20, 15))
        axes = axes.flatten()

        # Select a subset of systems to visualize (e.g., different nlin values)
        selected_systems = [s for s in systems if int(s.split()[4]) in [64, 96, 128, 192]]

        for idx, system in enumerate(selected_systems[:20]):
            sets = np.array(system.split(" ")).astype(int)
            L = sets[0]
            mring = sets[1]
            nring = sets[2]
            mlin = sets[3]
            nlin = sets[4]
            
            # Run a single trial to get the final graph
            res, G = run_trial(mring, nring, mlin, nlin, max_steps=args.max_steps, seed=42)
            
            # Calculate degree distribution
            degrees = [G.degree(n) for n in G.nodes()]
            max_degree = max(degrees) if degrees else 0
            
            # Count fraction of nodes with each degree
            degree_counts = np.bincount(degrees, minlength=int(max_degree) + 1)
            degree_fractions = degree_counts / len(G.nodes()) if G.number_of_nodes() > 0 else degree_counts
            
            k_values = np.arange(len(degree_fractions))
            
            axes[idx].bar(k_values, degree_fractions, edgecolor='k', alpha=0.7)
            axes[idx].set_xlabel('Degree (k)')
            axes[idx].set_ylabel('Fraction of nodes')
            axes[idx].set_title(f'nring={nring}, nlin={nlin}')
            axes[idx].grid(alpha=0.3)

        plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'degree_distributions.png'), dpi=150)
        plt.show()
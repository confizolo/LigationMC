# DNA ligation

# For the topology_reconstruction:
Python codes for the topology reconstruction after polymer condensation. The codes process files output from LAMMPS only. Namely, a bond list bond.* is read as input. Example input LAMMPS script (run.lmp) is also provided in this repository. 

==================== input and output files =======================================

input files: bond.* lammps bond lists.

output files: hist_linear_*.txt, hist_ring_*.txt: histograms of the polymer lengths in linear and ring state vs time. 

average_length.txt: the average polymer length in time (average over linear and rings)

Stat_polys.txt: first column = time, second column = the number of polymers in linear state, third column = the number of polymers in ring state

===================== technical details ============================================

To run: with python 3 interpreter run the main.py file. The input parameters are Nmols = initial number of molecules in the system before the condensation process begins, Nmon = initial number of beads per polymer, t_0 = the first time frame to read (i.e., bonds.0 if t_0 = 0), t_max = the last time frame to read (i.e., bonds.1000 if t_max = 1000), dt = the interval between two time frames (how often the bond lists are saved from LAMMPS), sim_dt = the timestep used during the simulation. 

==================== information on the codes ======================================

main.py: The main script that reads the input parameters and calls the readoutput_lmp.topo_reconstruct.

readoutput_lmp.py: Topology reconstruction code. The bond list are read here and the lookup subrutine is called to find the connectivity matrix between all the connected beads. 

properties.py: some properties, such as the polydispersity index, are calculated here. More properties can be added in this script later on by the user. 

# For the generalised_smoluchowski:

MATLAB codes for the numerical solution of the Smoluchowski polymer condensation equations with an additional sink term to accound for ring formation. The explicit Euler scheme is used for the integration of the population equations. 

==================== input and output files =======================================

input files: average_length.txt. This file is generated after the topology reconstruction (see topology_reconstruction above).

output files: none

output variables: topological number kappa, a dimensionless number related to the ratio of ring rate formation/ linear rate formation. The ring and linear rates account for the linear and ring polymer chains forming during a condesation process. 

===================== technical details ============================================

To run: In MATLAB version 2020a (or 2019) run the smol_fitting.m file. The input parameters are Nframes = number of lines in the average_length.txt files vol = volume used during the simulation (with LAMMPS or any other software). 

==================== information on the codes ======================================

smol_fitting.m: The main script where the least squares fitting (LSQ) is performed on the average_length calculated after a simulation where polymer condensation was happening. The topological parameter kappa is exported.

Obj_smoluchowski.m: The objective function used for the least squares fitting. This function is custom made and calls the exEuler_smoluchowski.m in every LSQ iteration. 

exEuler_smoluchowski.m: Explicit Euler scheme for solving the ordinary differential equations of the generalised Smoluchowski equation. 

# To run simulations go to the LAMMPS directory

run2.lam is an *example* to run MD simulations undergoing stochastic ligation of polymer ends 

melt... is an *example* initial configuration

reconstruct_topology is a mathematica script to reconstruct the topology of the system 

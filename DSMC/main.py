#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Dr Maria Panoukidou
Affiliation: The University of Edinburgh
"""

import readoutput_lmp
import pathlib
import numpy as np
import os;
import subprocess 

#### parameters ####
Nmols=200;  #number of molecules in the simulation @ the initial set up
Nmon=174; #number of beads per chain in the simulation @ the initial set up
t_0 = 0; #first frame to read
t_max = 1000000; #last frame to read
dt = 1000; #interval between frames
sim_dt = 0.01; #timestep used in simulation for the integration of equations of motion
h = (t_max-t_0)/dt + 1;
#####################

path=pathlib.Path().absolute()   #take current path

for i in range(0,40):  #loop over replicas

    #the lammps output files to read are located in folders named data_i where i represents a replica. 
    path2 = str(path)+'/data'+str(i+1)
    # delete any pre-existing folder with the previous results
    args = ["rm -rf",'"'+str(path2)+"/output"+'"']
    cmd = " ".join(args)
    subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True) == 0
    #create a folder for the output files after the topology reconstruction
    args = ["mkdir",'"'+str(path2)+"/output"+'"']
    cmd = " ".join(args)
    subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True) == 0
    path_out = os.path.join(path2, "output");
    Chains=[];
    Rings=[];
    
    #call the module that reads the bond lists to reconstruct topology
    readoutput_lmp.topo_reconstruct(Chains,Rings,Nmols,Nmon,t_0,t_max,dt,sim_dt,str(path2)+"/",str(path_out)+"/")
    num_polys = [row[0] for row in Chains];
    num_rings= [row[1] for row in Rings];
    tt = np.arange(h+1)
    #store data in the Stat_polys.txt file with info on the number of rings and linear chains vs time
    with open(os.path.join(path_out,"Stat_polys.txt"), "a") as myfile:
        for i in range(len(num_polys)):
            myfile.write("%f  %f  %f\n" % (tt[i], num_polys[i],num_rings[i]))

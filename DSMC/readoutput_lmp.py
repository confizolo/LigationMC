#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Dr Maria Panoukidou
Affiliation: The University of Edinburgh
"""

import os;
import properties as prop;

def topo_reconstruct(Chain_Data,ring_data,Nmols,Nmon,t0,tmax,dt,sim_dt,path,pathout):

     #open files to read
     timestep = t0;
     index = 1; 
            ##(1) Read bond list per timestep
                                   
            #lists of variables in time
            #Chain_Data=[]; #[number of linear chains]
            #ring_data=[]; #[Rings list, num of ring chains, minimum length of ring, maximum length of ring]
           
            #statistics lists
            #Stat_NumberBeads=[]; #number of beads per poly
            #Stat_BeadsInRings=[]; #number of beads per ring

     while timestep <=tmax:
        if(timestep == 0): 
             fullfilename = path+'bonds.0';
        else: 
             fullfilename = path+'bonds.'+str(timestep)+'000';
         
        with open(fullfilename, 'r') as bondFile:
            #read n (number of bonds in total, after the ligation)
            for i in range(9):
                line = bondFile.readline();
                if(i==3): n = int(line);
            
            #read bond connections
            Con = [];
            for i in range(n):
               line = bondFile.readline().split();           
               Con.append([int(line[0]),int(line[1]),int(line[2])]);
            
            ##find rings and linear chains
            Polys=[]; Rings =[];  
            PolyRings=[];
            m1=[]; m2=[];
            for i in range(len(Con)):         
                m1.append(Con[i][1])
                m2.append(Con[i][2])
            # the lookup module searches for new connections between beads and reconstructs the topology
            # Polys list contains the lists of beads indexes that are connected in linear polymers, either pre-existing or formed after the ligation
            # PolyRings list contains the lists of beads indexes that are connected in ring polymers formed after the ligation
            if len(Con)>0: Rings = lookup(Polys,PolyRings,m1,m2);
            
            # at the final time frame only, read the lists of linear and ring polymers and calculate the polydispersity index of the system
            # the index is then stored in the PDI.txt file
            if timestep==tmax:
                pdi = prop.PDI(Polys,PolyRings);
                with open(os.path.join(pathout,"PDI.txt"), "a") as myfile0:
                    myfile0.write("%f\n" % (pdi))
            
            #analysis of chains
            Stat_NumberBeads=[];
            for poly in Polys:
                # save the length of each linear polymer in the Stat_NumberBeads list. Length measured in beads.
                Stat_NumberBeads.append(len(poly));
            # save the number of linear chains, if there are any, in the Chain_Data list
            if len(Polys)>0:
                Chain_Data.append([len(Polys)]);
            else:
                Chain_Data.append([0]);
            # save the Rings list, num of ring chains, minimum length of ring, maximum length of ring in the ring_data list
            if len(PolyRings)>0:
                ring_data.append([Rings,len(PolyRings),min(map(len, PolyRings)),max(map(len, PolyRings))]);
            else:
                ring_data.append([0,0,0,0]);

            
            #analysis of rings
            Stat_BeadsInRings=[];
            for poly in PolyRings:               
                # save the length of each ring polymer in the Stat_BeadsInRings list. Length measured in beads.
                Stat_BeadsInRings.append(len(poly));   
            print('timestep '+str(timestep)+' done!!!!!!!!!');  


#################### HISTOGRAMS + SAVE IN FILES#########################
            bins=[];bins2=[];
            total_mols=0;
            ave_l=0;
            if len(Stat_BeadsInRings)>0: 
             for i in range(Nmon,max(Stat_BeadsInRings)+Nmon,Nmon): bins2.append(i);
             counts2 = [];
             for el in range(len(bins2)):
                 counts2.append(0);
             indx=0;
             for n in range(len(bins2)):
                matches = [x-1 for x in Stat_BeadsInRings if (x-1)==bins2[n]];
                counts2[indx] = len(matches);
                indx=indx+1;
             with open(pathout+'hist_rings_'+str(index)+'.txt','w') as myfile:
                 for oo in range(len(counts2)):
                     myfile.write("%f  %f \n" % (bins2[oo], counts2[oo]))
                     total_mols = total_mols + counts2[oo];
            if len(Stat_NumberBeads)>0:
             for i in range(Nmon,max(Stat_NumberBeads)+Nmon,Nmon): bins.append(i);
             counts = [];
             for el in range(len(bins)):
                 counts.append(0);  
             indx=0;
             for n in range(len(bins)):
                 matches = [x for x in Stat_NumberBeads if x==bins[n]];
                 counts[indx] = len(matches);
                 indx=indx+1;
             with open(pathout+'hist_linear_'+str(index)+'.txt','w') as myfile:
                  for oo in range(len(counts)):
                      myfile.write("%f  %f \n" % (bins[oo], counts[oo]))
                      total_mols = total_mols+ counts[oo];
             ave_l = (Nmols*Nmon)/total_mols; 
             with open(pathout+'average_length.txt','a') as myfile:
                  myfile.write("%.1f  %.10f \n" % (timestep*sim_dt, ave_l))
         
            timestep=timestep+dt; #udpate time for next iteration
            index = index+1;

            

########## LOOKUP FOR LINEAR POLYS OR RINGS #################


def lookup(Polys,PolyRings,m1,m2):
    #reformulate the ids
    modIds = list(set(m1+m2));
    
    # Create a dictionary for the indexes
    DicIds={};
    for i in range(len(modIds)): DicIds[modIds[i]]=i; 
      
    #create 0 based m vectors
    m1f=[];
    m2f=[];
    for i in range(len(m1)):
        m1f.append(DicIds[m1[i]]);
        m2f.append(DicIds[m2[i]]);
  
    #print('m1f m2f is:');
    #print(m1f);
    #print(m2f);
    

    #create connectivity matrix c where the rows are theparticle index,
    # the first column has the index of the left connected particle and the second column has the 
    # index of the right connected particle. The last column if a flag (0 or 1) and checks for the particles that 
    # are already read and belong to a molecule. 
    c=[];
    for i in range(len(modIds)): c.append([-1,-1,0])
    for i in range(len(m1f)):
        if c[m1f[i]][0]==-1:         c[m1f[i]][0]=m2f[i];
        elif c[m1f[i]][1]!=m2f[i]:  c[m1f[i]][1]=m2f[i];

        if c[m2f[i]][0]==-1: c[m2f[i]][0]=m1f[i];
        elif c[m2f[i]][0]!=m1f[i]: c[m2f[i]][1]=m1f[i];

    #print('c is:');
    #print(c);

    #create linear chain polymers by reading the connectivity matrix c
    chain=[];
    for i in range(len(modIds)):
        if c[i][2]==0 and c[i][1]==-1:
            chain.append(i);
            chain.append(c[i][0]);
            c[i][2]=1;
            previousMol=i;
            currentMol=c[i][0];
            c[currentMol][2]=1;
            while c[currentMol][1]!=-1:
                if c[currentMol][0]==previousMol:
                    chain.append(c[currentMol][1]);
                    previousMol=currentMol;
                    currentMol=c[currentMol][1];          
                    c[currentMol][2]=1;
                else:
                    chain.append(c[currentMol][0]);
                    previousMol=currentMol;
                    currentMol=c[currentMol][0];          
                    c[currentMol][2]=1; 
              
            Polys.append(chain);     
            chain=[];
    
    # Similarly to above, create Ring chains by looping over the matrix c
    chain=[];
    NmolsInPolys=0;
    for i in range(len(Polys)): 
         NmolsInPolys= NmolsInPolys+len(Polys[i]);
    
    if NmolsInPolys<len(c):
        for i in range(len(c)):
            if c[i][2]==0:
                RingStart=i;
                chain.append(i);
                chain.append(c[i][0]);
                c[i][2]=1;
                previousMol=i;
                currentMol=c[i][0];
                c[currentMol][2]=1;
                while currentMol!=RingStart:
                    if c[currentMol][0]==previousMol:
                        chain.append(c[currentMol][1]);
                        previousMol=currentMol;
                        currentMol=c[currentMol][1];
                        c[currentMol][2]=1;
                    else:
                        chain.append(c[currentMol][0]);
                        previousMol=currentMol;
                        currentMol=c[currentMol][0];
                        c[currentMol][2]=1;
                PolyRings.append(chain);
                chain=[];
         
    Rings=len(c)-NmolsInPolys;
        
    #remap the particle indexes to the initially assigned indexes by lammps.     
    for i in range(len(Polys)):  
        for z in range(len(Polys[i])):    
            Polys[i][z]=modIds[Polys[i][z]]; 
    for i in range(len(PolyRings)):
        for z in range(len(PolyRings[i])):
            PolyRings[i][z]=modIds[PolyRings[i][z]]
    return Rings;

# -*- coding: utf-8 -*-
"""
Author: Dr Maria Panoukidou
Affiliation: The University of Edinburgh
"""

import numpy as np;

def PDI(Polys,PolyRings):
    lens = [];
    lensR = [];
    for ring in PolyRings:
        lensR.append(len(ring)-1);
    
    for poly in Polys:
        lens.append(len(poly));
    
    n = [];
    N = [];
    X = np.array(lens);
    Y = np.array(lensR);
    for i in range(min(lens),max(lens)+1):
        indx = np.where(X==i)[0];
        if len(indx)>0:
            n.append(len(indx));
            N.append(i);
            
    for i in range(min(lensR),max(lensR)+1):
        indx = np.where(Y==i)[0];
        if len(indx)>0:
            n.append(len(indx));
            N.append(i);


    sqN = [N ** 2 for N in N];
    Mn = np.dot(n,N)/sum(n);
    Mw = np.dot(n,sqN)/np.dot(n,N);
    pdi = Mw/Mn;

    return pdi

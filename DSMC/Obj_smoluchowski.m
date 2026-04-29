%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Objective function for the lsqcurvefit that calculates       %
%            the average length of linear and ring condensated polymers   %
%            according to the Smoluchowski coagulation formula.           %
%                                                                         %
%                    Author: Dr Maria Panoukidou                          %
%                    The University of Edinburgh                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ydata = Obj_smoluchowski(k,xdata,vol)

% the array k contains the interaction coefficients k1 and k2 for linear
% and ring chain formation
k1 = k(1);
k2 = k(2);

alpha = 1;
nu = 0.6;
% array L with the lengths of polymers that can be formed/found in the
% system. In this particular case we have N = 174 monomers per polymer and
% 200 molecules, thus the longest polymer that can be formed will have
% length 34800 beads. 
L=174:174:34800;

% in this case we take the interaction rates as constant, however one could
% also assing different rate k1, k2 according to the length of the polymer.
% In our case this is taken into account in the exEuler_smoluchowski
% function, where the DeGenes theory is considered. 
K1=ones(length(L),length(L)).*k1;
K2=ones(length(L),1).*k2; 

% initial polymer number densities 
% n_t0_l is the number density of linear molecules with length l.
% n_t0_r is the number density of ring molecules with length l.
n_t0_l=[200;zeros(length(L)-1,1)]./vol;
n_t0_r=zeros(length(L),1)./vol;

% timestep dt
dt=xdata(2) - xdata(1);
% the total number of frames run in the simulation
Tsteps = length(xdata);

n_tp_L=n_t0_l;
n_tp_R=n_t0_r;
for Nstep=1:Tsteps-1
    % call exEuler_smoluchowski function which calculates the new molecule
    % number densities based on the rates k1 and k2
   [n_tn_L,n_tn_R]=exEuler_smoluchowski(n_tp_L,n_tp_R,dt,K1,K2,L,alpha,nu);
   
   % calculate the average linear chains' length at the current time frame
   % Nstep
   Lav_linear(Nstep)=sum(n_tn_L.*L')/sum(n_tn_L);  
   
   % calculate the average ring chains' length at the current time frame
   % Nstep
   if sum(n_tn_R)==0 
        Lav_ring(Nstep)=0;
   else
        Lav_ring(Nstep)=sum(n_tn_R.*L')/sum(n_tn_R); 
   end
   
   % calculate the total average length at the current time frame Nstep
   Lav_Total(Nstep)=(sum(n_tn_L.*L')+sum(n_tn_R.*L'))/(sum(n_tn_R)+sum(n_tn_L));
   
   % assing the new number densities to the n_tp variables and give them as
   % input to exEuler_smoluchowski in the next step.
   n_tp_L=n_tn_L; 
   n_tp_R=n_tn_R; 
end

% the final average length arrays are stored and the total average length
% is assigned to the output array ydata. If there are no rings formed in
% the system then the ydata will be simply the average linear chain lenght.
Lav_linear = [174,Lav_linear];
Lav_ring = [0,Lav_ring];
Lav_Total = [174,Lav_Total];

if(sum(Lav_ring~=0))
    ydata = Lav_Total;
else
    ydata = Lav_linear;
end
end
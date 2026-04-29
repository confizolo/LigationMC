%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Explicit Euler scheme for the numerical solution of          %
%            the Smoluchowski coagulation equation including              %
%            a ring formation term.                                       %
%                                                                         %
%            Author: Dr Maria Panoukidou                                  %
%            The University of Edinburgh                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n_t1_l,n_t1_r] = exEuler_smoluchowski(n_t0_l,n_t0_r,dt,K1,K2,L,alpha,nu)
n_t1_r=zeros(length(n_t0_r),1);  % ring chains array initialisation
n_t1_l=zeros(length(n_t0_l),1);  % linear chains array initialisation


%intermediate terms of smoluchowski equation
for k=1:length(L) 
  n_t1_l(k)=n_t0_l(k); 
  n_t1_r(k)=n_t0_r(k);
      for i=1:length(L)
          % 2nd term for molecules lost due to ligation bcause they formed bigger molecules
          if i+k<=length(L)
                n_t1_l(k)=n_t1_l(k)-dt*n_t0_l(i)*n_t0_l(k)*K1(i,k)*(i^(-alpha)+k^(-alpha))*(i^(nu)+k^(nu));
          end
          % 1st term of equation for molecules forming length k by addition
          for j=1:length(L)
              if i+j==k
                    n_t1_l(k)=n_t1_l(k)+dt*0.5*n_t0_l(i)*n_t0_l(j)*K1(i,j)*(i^(-alpha)+j^(-alpha))*(i^(nu)+j^(nu));
              end
          end       
      end   
    
    % 3rd term of the equation accounting for the ring formation.
    % The rings formed will be subtructed by the linear chains population.
    n_t1_r(k)=n_t1_r(k)+dt*K2(k)*k^(-4*nu)*n_t0_l(k);
    n_t1_l(k)=n_t1_l(k)-dt*K2(k)*k^(-4*nu)*n_t0_l(k);
      
end



end
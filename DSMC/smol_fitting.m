%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Matlab script for least squares fitting of the simulation    %
%            data according to the Smoluchowski coagulation formula.      %
%                                                                         %
%                    Author: Dr Maria Panoukidou                          %
%                    The University of Edinburgh                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%%  READ FILES AND TAKE AVERAGES OF LENGTH
tic
% take the current path
P1 = pwd;
% parameters
Nframes = 1001;
alpha = 1;
nu = 0.6;
vol = (346.938)^3; % WARNING: the volume should change according to the concentration!!!

% read files with data of the average length calculated after the topology
% reconstruction and take an average every 10 replicas. File format: time , length
filename = 'average_length.txt';
j = 1;
n = 0;
len = zeros(Nframes,1);
for i = 1:40 % loop over the replicas
    folder = ['data',num2str(i),'/output'];
    fullname = fullfile(P1,folder,filename);
    if isfile(fullname)
        f1 = load(fullname);
        if(size(f1,1)>=Nframes)
            len = len + f1(1:Nframes,2);
            n = n + 1;
        else
            continue;
        end
    end
    if(mod(i,10)==0)
        len = len./n;
        eval(sprintf('l_%d=len;',j));
        j = j + 1;
        len = zeros(Nframes,1);
        n = 0;
    end
end
t = f1(1:Nframes,1);


%% FITTING 

my_obj = @(k,xdata) Obj_smoluchowski(k,xdata,vol);
% initial guess for the reates, [K1:linear chain rate, K0:ring chain rate]
x0 = [10^(-6), 10^(-8)];
% lower and upper bounds for the rates
lb=[0; 0];
ub=[40; 2];
% least squares fitting
x = lsqcurvefit(my_obj,x0,t,l_1',lb,ub);
x1 = lsqcurvefit(my_obj,x0,t,l_2',lb,ub);
x2 = lsqcurvefit(my_obj,x0,t,l_3',lb,ub);   
x3 = lsqcurvefit(my_obj,x0,t,l_4',lb,ub);

% the rates found by the fitting are asigned to the matrix rates
rates(1,:) = x;
rates(2,:) = x1;
rates(3,:) = x2;
rates(4,:) = x3;

% mean values and standard deviations are calculated
m_k1(1,1) = mean(rates(:,1));
m_k1(1,2) = std(rates(:,1));
m_k0(1,1) = mean(rates(:,2));
m_k0(1,2) = std(rates(:,2));
ratio_std = std(2.*rates(:,2)./((200/vol).*rates(:,1)));

% kappa = 2*k_o/(n*k_1) where n = molecules/vol = 200/vol
Y = m_k0(1,1)./m_k1(1,1);
A = 2./(200./vol);
kappa = Y.*A';

toc
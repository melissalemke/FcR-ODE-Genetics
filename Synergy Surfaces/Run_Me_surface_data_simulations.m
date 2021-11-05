% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate -- Edited by Robby Theisen, 3D Isobologram
% Last edit: August 9, 2021

% Performs the simulations needed to make a combined and individually
% changed surface. Outputs:

%yfull = 50x50x33 matrix 1st & 2nd dimension are the changes in p1 and p2, 
            %3rd is complex formation at each of the p1 and p2 combinations
%ybase = vector of compelex concentrations at steady state
%y_p1 & y_p2 = vector of complex 33 concentration at steady state with
            %change in p1 or p2 at each of the dParams multipliers 
%params = baseline parameters
%p1 & p2 = the parameter indeces of the two params you are changing 
%dParams_p1 &dParams_p2 = the parameter modifiers that will be multiplied 
            %by the params at p1 and p2 
%complexes = name of the 33 output complexes recorded
%paramnames = name of the 22 input parameters
%fcr_names = names of the FcRs in order

tic
clear;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
fcr_id = 1;

% Read in parameters and ICs
file_to_load = Parameters(fcr_id);
load(file_to_load)

dParams_p1 = logspace(log10(0.1),log10(40),50); %Range for affinity
dParams_p2 = logspace(log10(.02),log10(3),50);  %Range for concentration
%makes vector of 50
% logarithmically spaced values between 0.02X and 45X to 
% multiply our baseline parameters by

% Run unaltered simulation
[ybase, steadystate] = Simulate(params, paramnames, [],[]);
if steadystate>0
    disp(["Baseline did not run to steady state"])
else
    disp(["Baseline ran to steady state"])
end
save(strcat('8_9_2021','_Base_Simulation_',fcr_ttl,'.mat'),'ybase','params',...
    'paramnames','complexes');

% Run full 2D sensitivity analysis
delete(gcp('nocreate')) %shutdown any parallel pools open
p1 = 9;%kon-IgG1-FcR index
p2 = 17;%IgG1 conc index
[yfull, filettl] = Sensitivity_2D(params, dParams_p1, dParams_p2,...
    paramnames,length(ybase(1,:)),p1,p2);

yfull = reshape(yfull, [length(dParams_p1),length(dParams_p1),length(yfull(1,:))]);

% Run 1D Sensitivity Analyses
y_p1 = Sensitivity_1D(params, paramnames, dParams_p1, p1);
y_p2 = Sensitivity_1D(params, paramnames, dParams_p2, p2);

save([erase(filettl,".mat")+"_"+fcr_ttl+".mat"])

toc
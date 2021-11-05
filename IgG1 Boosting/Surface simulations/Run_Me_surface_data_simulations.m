% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: October 11th, 2021

% Performs the simulations needed to make a surface. Outputs:

%yfull = 50x50x33 matrix 1st & 2nd dimension are the changes in p1 and p2, 
            %3rd is complex formation at each of the p1 and p2 combinations
%params = baseline parameters
%p1 & p2 = the parameter indeces of the two params you are changing 
%dParams = the parameter modifiers that will be multiplied by the params at
            %p1 and p2 
%complexes = name of the 33 output complexes recorded
%paramnames = name of the 22 input parameters
%fcr_names = names of the FcRs in order

tic
clear;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
% fcr_id = 1;
%load the data for FcR #1 to get FcR names
file_to_load = Parameters(1); 
load(file_to_load)

for fcr_id = 1:length(fcr_names)
    fcr_ttl = fcr_names(fcr_id);
% Read in parameters and ICs
file_to_load = Parameters(fcr_id); 
load(file_to_load)

dParams = logspace(log10(0.001),log10(20),50); %makes vector of 50
    % logarithmically spaced values between 0.004X and 20X to 
    % multiply our baseline parameters by

% Run unaltered simulation
        timecourse = false; % do not plot timecourse
% Run Baseline simulation
[ybase steadystate complexes] = Simulate(params, paramnames, complexes, fcr_ttl);

% Display whether or not simulation reached steady state as defined within
% "simulate" function
if steadystate>0
              disp(["Baseline did not run to steady state"])
        else
            disp(["Baseline ran to steady state"])
end
save(strcat(datestr(today()),'_Base_Simulation_',fcr_ttl,'.mat'),'ybase','params',...
    'paramnames','complexes');

% Run full 2D sensitivity analysis
delete(gcp('nocreate')) %shutdown any parallel pools open
p1 = 17;%IgG1 conc index
p2 = 19;%IgG3 conc index
[yfull, filettl] = Sensitivity_2D(params, dParams,...
    paramnames,length(ybase(1,:)),p1,p2,"",fcr_ttl);

yfull = reshape(yfull, [length(dParams),length(dParams),length(yfull(1,:))]);

save(filettl)

end

toc
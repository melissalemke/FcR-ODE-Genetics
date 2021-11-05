% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: October 11th, 2021

%% disclaimer: this simulation may take around 30 hours to complete
    % shorten the boosts vector or the number of rand_trials to run a shorter test

% Simulate IgG1 concentration and affinity to FcR boosts in all
% vaccinees, projected into each mixed allotype population. Mixed allotype
% populations are formed by randomly assigning the 105 vaccinees to an
% allotype based on the population makeup of each population, and
% converting thier initial IgG concentration accordingly. We assume all
% vaccinees start as G1m1,3. This random assignment is perform 25 times to
% form a more robust pooled population of the given population makeup. So
% for each pooled population n = 2,625

% The outputs to save in a mat file are the following, 
% each population saved into its own .mat file:

%all_run_rand = 1 cells for each random trial - each contain: 2x11x4x105x33 matrix
    %1st dimension is boost type: concentration or affinity
    %2nd is boost level
    %3rd is FcR
    %4th is person
    %5th is complexes
%param_idv_rand = 1 cells for each random trial - each contain: 2x11x4x105x22 matrix
    %1st dimension is boost type: concentration or affinity
    %2nd is boost level
    %3rd is FcR
    %4th is person
    %5th is input parameters
%complexname = name of the 33 output complexes recorded
%paramnames = name of the 22 input parameters
%IgG_FcR_data = 105x7 matrix: 1st dimension is person 2nd is the personal
    %input concentrations in col 1-4 for IgG1-4, ignore col 5-7
%patient_id = person identifiers
%FcR_names = names of the FcRs in order
%pop_makeup = 10x3 matrix mix of allotypes in each population, 1st
    %dimension is population, 2nd is allotype
%allotypes = names of allotypes in order
%boost_names = names of the boosting types in order
%boosts = boosting levels vector
%og_allotype_ind = cell for each random trial containing the allotype 
    %indece for the original allotype assumed of each person 
%new_allotype_ind = cell for each random trial containing the allotype 
    %indece assigned to (and converted to) for each person in the trial
    
function Run_Me_personal_simulations()
clear();
tic()
%% choose which strain to simulate, A244 or BAL
strain = 'A244';

boost_names = ["k_{on} IgG1-FcR","IgG1 conc"];
boosts = [0 0.1 0.25 0.5 0.75 1 2 2.5 5 7.5 10];

file_name = ['RV144_',strain,'_Personal_Input_v1-v105.xlsx'];

%% allotype projection: RV144 data is from Thai people - we expect to have
% g1m1 and g1m1,3, so we are assuming half of the samples come from g1m1
% and half came from g1m1,3 - we will randomly assign samples to one of
% these and convert to the given allotype from there
%IgG1 IgG2 IgG3 IgG4
g1m1_conv = [0.43 3.0 1.8 0.22; %conversion from g1m1 to g1m-1,3
    1   1   1   1;  %conversion from g1m1 to g1m1
    1.4 2.4 6.7 0.011]; %conversion from g1m1 to g1m1,3

%IgG1 IgG2 IgG3 IgG4
g1m13_conv = [0.30 1.3 0.26 20; %conversion from g1m1,3 to g1m-1,3
    0.70 0.41 0.15 91; %conversion from g1m1,3 to g1m1
    1    1   1   1];  %conversion from g1m1,3 to g1m1,3

%IgG1 IgG2 IgG3 IgG4
g1m3_conv = [1    1   1   1; %conversion from g1m-1,3 to g1m-1,3
    2.3 0.33 0.57 4.5; %conversion from g1m-1,3 to g1m1
    3.4 0.80 3.8 0.049];  %conversion from g1m-1,3 to g1m1,3

conversions = {g1m3_conv, g1m1_conv, g1m13_conv};
allotypes = ["G1m-1,3" "G1m1" "G1m1,3"];

% populations = {'Asian';'Caucasian';'African'};
populations = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J'};
og_pop_id = 8;

%G1m-1,3 G1m1 G1m1,3
pop_makeup = [1 0 0; %A
    0.66 0.17 0.17; %B
    0.5 0 0.5; %C
    0.5 0.5 0;
    0.33 0.33 0.33;
    0.17 0.66 0.17;
    0.17 0.17 0.66;
    0 0 1;
    0 1 0;
    0 0.5 0.5];

% Read in personal data
indiv_data = readtable(file_name);

%Descibes column index of each desired parameter
IgG_FcR_data = indiv_data{:,2:end};

patient_id = indiv_data{:,1};

IgG_FcR_data(isnan(IgG_FcR_data))=0;%turn nans to zero in data
g1_idx = 1;
g2_idx = 2;
g3_idx = 3;
g4_idx = 4;

complexname=['env-IgG1', 'env-IgG2', 'env-IgG3', 'env-IgG4',...
    'env-IgG1-IgG1', 'env-IgG1-IgG2','env-IgG1-IgG3', 'env-IgG1-IgG4',...
    'env-IgG2-IgG2', 'env-IgG2-IgG3','env-IgG2-IgG4','env-IgG3-IgG3',...
    'env-IgG3-IgG4','env-IgG4-IgG4','FcR-env-IgG1-IgG1',...
    'FcR-env-IgG1-IgG2','FcR-env-IgG1-IgG3', 'FcR-env-IgG1-IgG4',...
    'FcR-env-IgG2-IgG2', 'FcR-env-IgG2-IgG3','FcR-env-IgG2-IgG4',...
    'FcR-env-IgG3-IgG3','FcR-env-IgG3-IgG4','FcR-env-IgG4-IgG4',...
    "IgG1", "IgG2", "IgG3", "IgG4", "env", "FcR",...
    "FcR complexes with IgG1 and IgG3 only", ...
    "FcR complexes including IgG1 and IgG3","All FcR complexes"];

% Different FcR affinity values(nM-1s-1), rows are IgG1-IgG4 Fc affinity,
% columns are the FcR (FcR3aV FcR3aF FcR2aH FcR2aR  respectively)
FcR_kon = [20	11.7	52  35; % IgG1-Fc kon
    0.7	0.3	4.5     1; % IgG2-Fc kon
    98	77	8.9     9.1; % IgG3-Fc kon
    2.5	2	1.7     2.1]*10^-6; % IgG4-Fc kon
FcR_names = ["FcR3a-V158" "FcR3a-F158" "FcR2a-H131" "FcR2a-R131"];

%% Run baseline simulations
clear all_run
clear param_idv

fcr_id = 0;
for fcr_idx = 1:length(FcR_kon(1,:)) % columns in the spreadsheet
    fcr_id = fcr_id + 1;
    for p_num = 1:length(patient_id)
        
        % get personal paramters for the given FcR
        [params,paramnames,~] = ...
            Parameters_indiv_FcR(FcR_kon(1,fcr_id), FcR_kon(2, fcr_id), ...
            FcR_kon(3, fcr_id), FcR_kon(4, fcr_id), ...
            IgG_FcR_data(p_num,g1_idx), IgG_FcR_data(p_num,g2_idx), ...
            IgG_FcR_data(p_num,g3_idx), IgG_FcR_data(p_num,g4_idx),boost_names(1),...
            1);
        
        % run the simulation
        [ybase, steadystate] = Simulate(params,paramnames,[],[]);
        
        % display if steady state was reached or not (defined in
        % "simulate.m"
        if steadystate>0
            disp(["Patient "+p_num+" did not reach steady state"])
        else
            disp(["Patient "+p_num+" at steady state"])
        end
        % save data in the larger matrix
        all_run_rv144(fcr_id, p_num, :) = ybase;
        param_idv_rv144(fcr_id, p_num,:) = params;
    end
end
%% Run random simulations to mimic a different population
for pop_id = 1:length(populations)
    for rand_trial = 1:25
        %% randomly assigning allotypes based on the population makeup
        og_allotype_ind{rand_trial} = get_allotype(pop_makeup, og_pop_id,patient_id);
        new_allotype_ind{rand_trial} = get_allotype(pop_makeup, pop_id,patient_id);
        
        %% simulations
        for boost_type = 1:length(boost_names)
            for boost_level = 1:length(boosts)
                fcr_id = 0;
                for fcr_idx = 1:length(FcR_kon(1,:)) % columns in the spreadsheet
                    fcr_id = fcr_id + 1;
                    
                    for p_num = 1:length(patient_id)
                        %allotype conversion
                        og_allo = og_allotype_ind{rand_trial}(p_num);
                        new_allo = new_allotype_ind{rand_trial}(p_num);
                        cf = conversions{og_allo}(new_allo,:);
                        
                        % get personal paramters for the given FcR
                        [params,paramnames,~] = ...
                            Parameters_indiv_FcR((FcR_kon(1,fcr_id)),FcR_kon(2, fcr_id), FcR_kon(3, fcr_id), FcR_kon(4, fcr_id), ...
                            IgG_FcR_data(p_num,g1_idx)*cf(1),...
                            IgG_FcR_data(p_num,g2_idx)*cf(2), ...
                            IgG_FcR_data(p_num,g3_idx)*cf(3), ...
                            IgG_FcR_data(p_num,g4_idx)*cf(4),...
                            boost_names(boost_type),boosts(boost_level));
                        
                        % run the simulation
                        [ybase, steadystate] = Simulate(params,paramnames,[],[]);
                        
                        % display if steady state was reached or not (defined in
                        % "simulate.m"
                        if steadystate>0
                            disp(["Patient "+p_num+" did not reach steady state"])
                        else
                            disp(["Patient "+p_num+" at steady state"])
                        end
                        % save data in the larger matrix
                        all_run(boost_type,boost_level,fcr_id, p_num, :) = ybase;
                        param_idv(boost_type,boost_level,fcr_id, p_num,:) = params;
                    end
                end
            end
        end
        all_run_rand{rand_trial} = all_run;
        param_idv_rand{rand_trial} = param_idv;
    end
    save([strain,'Population ',populations{pop_id},'_boost simulations_',datestr(today()),'.mat'],...
        'complexname','paramnames', 'IgG_FcR_data',...
        'patient_id','FcR_names','all_run_rand','param_idv_rand','allotypes',...
        'og_allotype_ind','new_allotype_ind','boost_names','boosts',...
        'populations','pop_makeup');
end
toc()
end

function allotype_ind = get_allotype(pop_makeup, pop_id,patient_id)
%num_allotypes = sum(pop_makeup(pop_id,:)); %how many allotypes in the original population
num_patients_each_allotype = pop_makeup(pop_id,:)*length(patient_id); %patients we assume have each allotype in the og population

%randomly round the # in each allotype so they are all integers
if sum(num_patients_each_allotype~=round(num_patients_each_allotype))
    pos_ind = find(num_patients_each_allotype);
    add_ind_in_pos_ind = randperm(length(pos_ind),1);
    ind_to_add = pos_ind(add_ind_in_pos_ind);
    
    num_patients_each_allotype = floor(num_patients_each_allotype);
    num_patients_each_allotype(ind_to_add) = num_patients_each_allotype(ind_to_add)+1;
    
    if sum(num_patients_each_allotype) ~= length(patient_id)
        unused_ind = find(pos_ind ~= add_ind_in_pos_ind);
        rand_unused_ind = randperm(2,1);
        next_ind_to_add = unused_ind(rand_unused_ind);
        num_patients_each_allotype(next_ind_to_add) = num_patients_each_allotype(next_ind_to_add)+1;
    end
    if sum(num_patients_each_allotype) ~= length(patient_id)
        last_unused_ind = find(unused_ind ~= next_ind_to_add);
        last_ind_to_add = unused_ind(last_unused_ind);
        num_patients_each_allotype(last_ind_to_add) = num_patients_each_allotype(last_ind_to_add)+1;
    end
end

% randomly assign patients to the allotypes present
allotype_ind = zeros(length(patient_id),1); %index, 1 = G1m-1,3, 2 = G1m1, 3 = G1m1,3

g1m3_ind = randperm(leng  th(patient_id),num_patients_each_allotype(1));
allotype_ind(g1m3_ind) = 1;
ind_left = find(allotype_ind==0); % this ensures each person gets an assignment

g1m1_ind = ind_left(randperm(length(ind_left),num_patients_each_allotype(2)));
allotype_ind(g1m1_ind) = 2;
ind_left = find(allotype_ind==0);

g1m13_ind = ind_left(randperm(length(ind_left),num_patients_each_allotype(3)));
allotype_ind(g1m13_ind) = 3;
end
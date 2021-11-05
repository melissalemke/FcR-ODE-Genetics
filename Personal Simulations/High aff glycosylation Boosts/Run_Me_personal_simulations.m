% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: October 11th, 2021

% Simulate baseline and a high affinity glycosylation boost in all 
% vaccinees, projected into each allotype: The outputs to save in a mat 
% file are the following:

%all_run_rand = 1 cells for each population - each contain: 4x2x105x33 matrix 
            %1st dimension is FcR, 2nd is boost level, 3rd is person, 
            %4th is complexes
%param_idv_rand = 1 cells for each population - each contain: 4x2x105x22 matrix 
            %1st dimension is FcR, 2nd is boost level, 3rd is person, 
            %4th is parameters
%complexname = name of the 33 output complexes recorded
%paramnames = name of the 22 input parameters
%IgG_FcR_data = 105x7 matrix: 1st dimension is person 2nd is the personal
            %input concentrations in col 1-4 for IgG1-4, ignore col 5-7
%patient_id = person identifiers
%FcR_names = names of the FcRs in order
%pop_makeup = 10x3 matrix mix of allotypes in each population, 1st
            %dimension is population, 2nd is allotype
%populations = names of populations
%allotypes = names of allotypes
%boost_labels = names of the boosting
%IgG1_boost = vector of multiplicative boosts of IgG1
%kon_boost = vector of multiplicative boosts of kon IgG1-FcR

function Run_Me_personal_simulations()
clear();
tic()
%% choose which strain to simulate
strain = 'A244';
%% boosting levels - boosts are multiplied by baseline parameters
IgG1_boost = [1 1];
kon_boost = [1 31];
boost_labels = {'No_boost', 'High affinity glycosylation IgG1-FcR_kon_boost'};

file_name = ['RV144_',strain,'_Personal_Input_v1-v105.xlsx'];

%% allotype projection: RV144 data is from Thai people - we expect to have
% g1m1,3 and convert to the given allotype from there
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

populations = {'G1m-1,3';'G1m1';'G1m1,3'};
og_pop_id = 2;

          %G1m-1,3 G1m1 G1m1,3
pop_makeup = [1 0 0; %A
              0 1 0; 
              0 0 1]; 

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
            IgG_FcR_data(p_num,g3_idx), IgG_FcR_data(p_num,g4_idx),IgG1_boost(1),...
            kon_boost(1));
        
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
        %% randomly assigning allotypes based on the population makeup
        og_allotype_ind{pop_id} = get_allotype(pop_makeup, og_pop_id,patient_id);
        new_allotype_ind{pop_id} = get_allotype(pop_makeup, pop_id,patient_id);
       
        %% simulations
        fcr_id = 0;
        for fcr_idx = 1:length(FcR_kon(1,:)) % columns in the spreadsheet
            fcr_id = fcr_id + 1;
            for boost_ind = 1:length(IgG1_boost)
            for p_num = 1:length(patient_id)
                %allotype conversion
                og_allo = og_allotype_ind{pop_id}(p_num);
                new_allo = new_allotype_ind{pop_id}(p_num);
                cf = conversions{og_allo}(new_allo,:);
                
                % get personal paramters for the given FcR
                [params,paramnames,~] = ...
                    Parameters_indiv_FcR((FcR_kon(1,fcr_id)),FcR_kon(2, fcr_id), FcR_kon(3, fcr_id), FcR_kon(4, fcr_id), ...
                    IgG_FcR_data(p_num,g1_idx)*cf(1),...
                    IgG_FcR_data(p_num,g2_idx)*cf(2), ...
                    IgG_FcR_data(p_num,g3_idx)*cf(3), IgG_FcR_data(p_num,g4_idx)*cf(4),IgG1_boost(boost_ind),kon_boost(boost_ind));
                
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
                all_run(fcr_id,boost_ind, p_num, :) = ybase;
                param_idv(fcr_id,boost_ind, p_num,:) = params;
            end
            end
        end
        all_run_rand{pop_id} = all_run;
        param_idv_rand{pop_id} = param_idv;
    end
     save([strain,'_boost simulations_',datestr(today()),'.mat'],...
        'complexname','paramnames', 'IgG_FcR_data','IgG1_boost','kon_boost',...
        'patient_id','FcR_names','all_run_rand','param_idv_rand',...
        'populations','pop_makeup','allotypes','boost_labels');
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
        
        g1m3_ind = randperm(length(patient_id),num_patients_each_allotype(1));
        allotype_ind(g1m3_ind) = 1;
        ind_left = find(allotype_ind==0); % this ensures each person gets an assignment
        
        g1m1_ind = ind_left(randperm(length(ind_left),num_patients_each_allotype(2)));
        allotype_ind(g1m1_ind) = 2;
        ind_left = find(allotype_ind==0);
        
        g1m13_ind = ind_left(randperm(length(ind_left),num_patients_each_allotype(3)));   
        allotype_ind(g1m13_ind) = 3;
    end
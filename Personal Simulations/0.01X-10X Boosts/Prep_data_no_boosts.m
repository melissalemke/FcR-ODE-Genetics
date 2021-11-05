clear

load('A244_personal_G1m13_to_all_populations_21-Jul-2021.mat')
fcr_id = 1;
com_id = 33;

pooled_trials = zeros(length(all_run_rand_combined{1})*length(patient_id),length(all_run_rand_combined));

for i = 1:length(all_run_rand_combined)
    trials_cell = all_run_rand_combined{i};
    trial_vector = [];
    for j = 1:length(all_run_rand_combined{1})
    trials(:,j,i) = trials_cell{j}(fcr_id,:,com_id);
    trial_vector = [trial_vector;trials_cell{j}(fcr_id,:,com_id)'];
    trial_med(j,i) = median(trials_cell{j}(fcr_id,:,com_id));
    end
    pooled_trials(:,i) = trial_vector; 
end


clear;
filedir = dir(['A244_boost simulations_*']);
load(filedir.name);

com_id = 33;
fcr_id = 1;

for i = 1:length(all_run_rand)
    baseline_data(:,i) = squeeze(all_run_rand{i}(fcr_id,1,:,com_id))';
    high_affinity_data(:,i) = squeeze(all_run_rand{i}(fcr_id,2,:,com_id))';
end

change_in_com_form = high_affinity_data-baseline_data;

dot_plot_data = change_in_com_form;
med_dot_data = median(dot_plot_data);
xlabels = allotypes;

figure()
x_val_dot_plot = rand(size(dot_plot_data))+(1:length(xlabels))-0.5;
scatter(x_val_dot_plot,dot_plot_data)
hold on
line([0.5 1.5;1.5 2.5;2.5 3.5]',[med_dot_data' med_dot_data']')
xticks(1:length(xlabels))
xticklabels(xlabels)
ylabel(["Change in complex"; "formation with high";"affinity glycosyltaion (nM)"])

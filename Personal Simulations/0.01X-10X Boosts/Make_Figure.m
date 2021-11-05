clear;
filedir = dir(['A244Population A_boost simulations_*']);
load(filedir.name);

com_id = 33;
fcr_id = 1;

for pop_id = 1:length(populations)
    file_name = dir(['A244Population ',populations{pop_id},'_boost simulations_*.mat']);
    load(file_name.name);
    
    for fcr_id = 1:length(FcR_names)
    all_run_pop = zeros(0);
    
    %pool the 25 trials - all_run is a cell for each populaiton containing
    %2x10x2625, 1st dim = boost type, 2nd dim = boost level, 3rd dim is person
    for trial_id = 1:length(all_run_rand)
        all_run_pop = cat(3,all_run_pop,squeeze(all_run_rand{trial_id}(:,:,fcr_id,:,com_id)));
    end
    all_run{pop_id,fcr_id} = all_run_pop;
    
    %get the change in complex formation from baseline
    for boost_type = 1:length(boost_names)
        for boost_level = 1:length(boosts)
            change_pop(boost_type,boost_level,:) = squeeze(all_run_pop(boost_type,boost_level,:))-...
                squeeze(all_run_pop(boost_type,1,:));
            
        end
    end
    
    change{pop_id,fcr_id} = change_pop;
    
    %get the medians in each population for each boosting case to make
    %the heatmaps
    for boost_type = 1:length(boost_names)
        change_med{boost_type,fcr_id}(:,pop_id) = median(change_pop(boost_type,:,:),3);
    end
    
   
    end
end
%get ratio of median change in com formation for first boosting type:second
for fcr_id = 1:length(FcR_names)
    change_med{3,fcr_id} = change_med{2,fcr_id}./change_med{1,fcr_id};
end
 boost_names(3) = "Change with IgG1 concentration/affinity to FcR";
 
% Plot FcR3aV heatmaps with all populations
fcr_id = 1;
for i = 1:length(boost_names)
    figure()
    heatmap(populations,string((boosts)*100)+"%",change_med{i,fcr_id})
    xlabel("Population")
    ylabel("Boost Level")
    title(FcR_names(fcr_id)+" "+boost_names(i))
end

% Plot allotype+FcR3 polymorphism genotype heatmaps with only 100% populations
populations_id = [8 9 1]; %these are the 100% 1 allotype populations

j=1;
for allotype = populations_id
    for fcr_id = 1:length(FcR_names(1:2))
        xlabels(j) = [allotypes(find(pop_makeup(allotype,:)))+" "+FcR_names(fcr_id)];
        for boost_type = 1:length(boost_names)
        heatmap_data{boost_type}(:,j) = change_med{boost_type,fcr_id}(:,allotype);
        end
        dot_plot_data(:,j) = squeeze(all_run{allotype,fcr_id}(1,1,1:105));
        j=j+1;
    end
end

figure()
x_val_dot_plot = rand(size(dot_plot_data))+(1:length(xlabels))-0.5;
scatter(x_val_dot_plot,dot_plot_data)
xticks(1:length(xlabels))
xticklabels(xlabels)
ylabel("Complex formation (nM)")

for i = 1:length(boost_names)
    figure()
    heatmap(xlabels,string((boosts)*100)+"%",heatmap_data{i})
    xlabel("Population")
    ylabel("Boost Level")
    title(boost_names(i))
end
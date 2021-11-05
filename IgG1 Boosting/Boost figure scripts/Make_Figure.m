% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: October 5th, 2021

%% make sure you move your results from personal simulations and 
%% surface simulation into this folder!
clear
%choose the output complex 
com_id = 33;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
% fcr_id = 1;

% load personal data just to set the size of x1 and x2 - will load again
% later
personal_data_path = 'A244_personal_baseline_all_fcrs_addition*.mat';
d = dir(personal_data_path);
file_name = d.name;
load(file_name);

group_names = ["Post Vaccination" "Post 145 nM IgG1 Boost"];
adds = [0 145];
color = [0 0 0; 0 0 1; 1 0 0; 1 0 1; 0 1 1]; % group color
fcr_names = ["FcR3a-V158" "FcR3a-F158"];

x1 = nan(length(patient_id),length(group_names)*length(fcr_names));
x2 = nan(length(patient_id),length(group_names)*length(fcr_names));

%choose the strain
strain = 'A244';

for fcr_id = 1:length(fcr_names)
%load the surface data
d = dir(['*_2D_IgG1_conc*'+fcr_names(fcr_id)+'*']);
file_name = d.name;
load(file_name);

%shorten fcr_names loaded from 2D file to only FcR3
fcr_names = ["FcR3a-V158" "FcR3a-F158"];% "FcR2a-H131" "FcR2a-R131"];

disp([paramnames(p1)+" "+paramnames(p2)])

x = dParams*params(p1);%nM
y = dParams*params(p2);%nM
z = yfull(:,:,com_id);%nM
z_fcrs(fcr_id,:,:) = z;

%% Load the personal baseline data!!! Move it to this folder
% try % will only run if you've moved the personal baseline data into this folder
d = dir(personal_data_path);
file_name = d.name;
load(file_name);
com_name = complexname(com_id);

all_run_full = all_run;
all_run = all_run{1};%baseline

% Prep data
FcR_model = squeeze(all_run(:,:,com_id))';

%plot the surface
surfacefig(p1,p2,FcR_model,fcr_id,patient_id, com_name,...
fcr_names,param_idv{1},x,y,z,paramnames,"", group_names,color);
end
figure(1)

%Add the lines and dots for each person
for fcr_id = 1:length(fcr_names)
for group_ind = 1:length(group_names)
    addition = adds(group_ind);
    [zs_pers(fcr_id,group_ind,:)] = getzpers(p1,p2,fcr_id,patient_id,param_idv{group_ind},x,y,squeeze(z_fcrs(fcr_id,:,:)));
end
end

for group_ind = 1:2%length(group_names)
    addition = adds(group_ind);
    param_idv_g  = param_idv{group_ind};
for pers = 1:length(patient_id)
    pers_color = color(group_ind,:);
line([param_idv_g(fcr_id,pers,p1)' param_idv_g(fcr_id,pers,p1)'],[param_idv_g(fcr_id,pers,p2)',...
    param_idv_g(fcr_id,pers,p2)'],[min(squeeze(zs_pers(:,group_ind,pers)))' max(squeeze(zs_pers(:,group_ind,pers)))'],'Color',pers_color,...
    'Marker','.','MarkerSize',20)
end
end

legend(fcr_names)

%plot the dot plot showing change between polymorphisms before and after
%boosting
for i = 1:length(all_run_full)
diff_btwn_polys(:,i) = all_run_full{i}(1,:,com_id)-all_run_full{i}(2,:,com_id);
end
dot_plot_data = diff_btwn_polys;
xlabels = ["Post vaccination" "Post IgG1 boost"];
med_dot_data = median(dot_plot_data);

figure()
x_val_dot_plot = rand(size(dot_plot_data))+(1:length(xlabels))-0.5;
scatter(x_val_dot_plot,dot_plot_data)
hold on
line([0.5 1.5;1.5 2.5]',[med_dot_data' med_dot_data']')
xticks(1:length(xlabels))
xticklabels(xlabels)
ylabel(["Change in complex"; "formation from FcR3aF";"to FcR3aV (nM)"])

%% Surface fig
function [zs_pers] = surfacefig(p1,p2,FcR_model,fcr_id,patient_id, com_name, FcR_names,param_idv,x,y,z,paramnames, groups, group_name,color)
figure(1)

surf_color = [1 1 1;0 0 0];

%plot surface - > currently everyother line to prevent crowding
surface(x(1:2:end),y(1:2:end),z(1:2:end,1:2:end),'LineStyle','-','FaceAlpha',0.5,'FaceColor',surf_color(fcr_id,:))%
ylabel([paramnames(p2)+" (nM)"])
xlabel([paramnames(p1)+" (nM)"])
xs = [min(x)*1 max(x)*1];
ys = [min(y)*1 max(y)*1];
zs = [min(min(z))*1 max(max(z))*1];
xlim(xs)
ylim(ys)

zlabel(["Complex", "Formation (nM)"])
set(gca, 'YScale', 'log','XScale','log')%,'Zscale','log'
set(gca,'FontSize',8)
hold on
end

function [zs_pers] = getzpers(p1,p2,fcr_id,patient_id,param_idv,x,y,z)
for pers = 1:length(patient_id)
[xplace xind] = min(abs(x-(param_idv(fcr_id,pers,p1))));
[yplace yind] = min(abs(y-param_idv(fcr_id,pers,p2)));
zs_pers(pers) = z(yind,xind);
end
end
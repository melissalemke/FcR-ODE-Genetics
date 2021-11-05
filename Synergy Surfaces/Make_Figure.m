% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate -- Edited by Robby Theisen
% Last edit: August 9, 2021

% plots three IgG1 concentration vs kon IgG1-FcR surfaces with either
% combinded changes, individual changes added together, or the difference
% between the two

clear
close all;
%choose the output complex 
com_id = 33;

%Choose which FcR: FcR3aV = 1 FcR3aF = 2 FcR2aH = 3 FcR2aR = 4
fcr_id = 1;

%choose the strain
strain = 'A244';
%load the surface data
d = dir(['*8_9_2021*IgG1_conc*']);
file_name = d.name;
load(file_name);

% Quick Check
disp([paramnames(p1)+" "+paramnames(p2)]);

% Setup Axes
x = dParams_p1*params(p1)*10^-6;%mM-1 to nM-1
y = dParams_p2*params(p2)*10^6;%mM to nM
z = yfull(:,:,com_id)*10^6;%mM to nM

%% Getting Patient Lines (all traces will be output in nM)
%%% Might be worth making this into a function that can handle making the
%%% traces on both the x and y axes

% Median Allotype Concentrations
g1m13_medconc = params(17);
g1m3_medconc = params(17) * .3;
g1m1_medconc = params(17) * .7;

% G1M1,3
newP = params;
trace1 = 10^6 * Sensitivity_1D(newP, paramnames, dParams_p1, p1);
add_trace1 = 10^6 * (ybase(33) + y_p1);

% G1M3 Median
newP(17) = g1m3_medconc;
trace2 = 10^6 * Sensitivity_1D(newP, paramnames, dParams_p1, p1);
ybase2 = Simulate(newP, paramnames, [], []);
add_trace2 = 10^6 * (ybase2(33) + y_p1);

% G1M1
newP(17) = g1m1_medconc;
trace3 = 10^6 * Sensitivity_1D(newP, paramnames, dParams_p1, p1);
ybase3 = Simulate(newP, paramnames, [], []);
add_trace3 = 10^6 * (ybase3(33) + y_p1);

%% Getting Glycosylation form lines
% Load Affinities
load("IgG1_glyc_kons.mat");
newP = params;

% Glycosylation form 18
newP(p1) = glyc_kon{"FcgRIIIA-158V", "G18"};
trace_g18 =  10^6 *Sensitivity_1D(newP, paramnames, dParams_p2, p2);
ybase_g18 = Simulate(newP, paramnames, [], []);
add_traceg18 = 10^6 *(ybase_g18(33) + y_p2);

% Glycosylation form 7
trace_g7 =  10^6 *Sensitivity_1D(params, paramnames, dParams_p2, p2);
ybase_g7 = Simulate(params, paramnames, [], []);
add_traceg7 = 10^6 *(ybase_g7(33) + y_p2);

load('redblue.mat')

%% Normal Surface fig
figure(1);
surf(x(1:2:end),y(1:2:end),z(1:2:end,1:2:end)-ybase(33)*10^6,'LineStyle','-');%
ylabel(["Initial IgG1";"Concentration (nM)"]);
xlabel([paramnames(p1);" (nM^{-1}s^{-1})"]);
hold on;
colormap(cm)
% Plot allotype traces
plot3(x, g1m13_medconc*ones(1, length(x))*10^6, trace1-ybase(33)*10^6, 'LineStyle','-', 'LineWidth',4,'color','w');
plot3(x, g1m3_medconc*ones(1, length(x))*10^6, trace2-ybase(33)*10^6, 'LineStyle','-', 'LineWidth',4,'color','k');
plot3(x, g1m1_medconc*ones(1, length(x))*10^6, trace3-ybase(33)*10^6, 'LineStyle','-', 'LineWidth',4,'color',[128 128 128]/256);
% Plot glycosylation forms
plot3(glyc_kon{"FcgRIIIA-158V", "G18"}*ones(1, length(x))*10^-6, y, trace_g18-ybase(33)*10^6, 'LineStyle','-', 'LineWidth',4)
plot3(glyc_kon{"FcgRIIIA-158V", "G2"}*ones(1, length(x))*10^-6, y, trace_g7-ybase(33)*10^6, 'LineStyle','-', 'LineWidth',4,'color',[0 114 188]/256)
legend(["", "G1M1,3", "G1M3", "G1M1", "High affinity glycosylation", "Baseline Fc\gammaRIIIa-V^{158}"]);
zlabel(["Complex", "Formation (nM)"]);
set(gca, 'YScale', 'log','XScale','log');%,'Zscale','log'
set(gca,'FontSize',8);
title(fcr_ttl);
view(-20,35)
hold off;
title(["Combined Changes in";"Affinity and Concentration"])

z_1d = (y_p1 + y_p2') * 10^6; % Convert mM to nM

% %% Normalized 1D Sensitivity Surface
% This surface is normalized such that we subtract out the effect of running
% it at baseline for either variable (see ppt for equation)
figure(2);
z_1d_norm = z_1d - (2 * ybase(33) * 10^6); % Subtracting out the value from a baseline simulation (twice to account for the effect being 2x)
surf(x(1:2:end), y(1:2:end), z_1d_norm(1:2:end,1:2:end));
set(gca, 'XScale', 'log', 'YScale', 'log');
title(fcr_ttl);
hold on;
colormap(cm)
% Plot Allotypes
plot3(x, g1m13_medconc*ones(1, length(x))*10^6, add_trace1- (2 * ybase(33) * 10^6), 'LineStyle','-', 'LineWidth',4,'color','w');
plot3(x, g1m3_medconc*ones(1, length(x))*10^6, add_trace2- (2 * ybase(33) * 10^6), 'LineStyle','-', 'LineWidth',4,'color','k');
plot3(x, g1m1_medconc*ones(1, length(x))*10^6, add_trace3- (2 * ybase(33) * 10^6), 'LineStyle','-', 'LineWidth',4,'color',[128 128 128]/256);
% Plot Glycosylaton patterns
plot3(glyc_kon{"FcgRIIIA-158V", "G18"}*ones(1, length(x))*10^-6, y, add_traceg18- (2 * ybase(33) * 10^6), 'LineStyle','-', 'LineWidth',4)
plot3(glyc_kon{"FcgRIIIA-158V", "G2"}*ones(1, length(x))*10^-6, y, add_traceg7- (2 * ybase(33) * 10^6), 'LineStyle','-', 'LineWidth',4,'color',[0 114 188]/256)
legend(["", "G1M1,3", "G1M3", "G1M1", "High affinity glycosylation", "Baseline Fc\gammaRIIIa-V^{158}"]);
set(gca, 'XScale', 'log', 'YScale', 'log');
ylabel(["Initial IgG1";"Concentration (nM)"]);
xlabel([paramnames(p1);" (nM^{-1}s^{-1})"]);
zlabel(["Complex", "Formation (nM)"]);
title(fcr_ttl);
title(["Affinity and Concentration","Changed Individually"])
view(-20,35)
z_lim(2) = max(max([z,z_1d_norm]));
z_lim(1) = min(min([z,z_1d_norm]));
zlim([-2 10])
z_norm_sub_norm = ((z - (ybase(33) * 10^6)) - z_1d_norm);

%% Adding Allotype Lines
figure(3);

surf(x(1:2:end), y(1:2:end), z_norm_sub_norm(1:2:end,1:2:end));
set(gca, 'XScale', 'log', 'YScale', 'log');
ylabel(["Initial IgG1";"Concentration (nM)"]);
xlabel([paramnames(p1);" (nM^{-1}s^{-1})"]);
zlabel(["Complex", "Formation (nM)"]);
title(fcr_ttl);
hold on;
colormap(cm)
% Calculate and plot the synergy 'slice' at each allotype concentration
syn_trace1 = (trace1 - ybase(33)*10^6) - (add_trace1 - 2*ybase(33)*10^6);
syn_trace2 = (trace2 - ybase(33)*10^6) - (add_trace2 - 2*ybase(33)*10^6);
syn_trace3 = (trace3 - ybase(33)*10^6) - (add_trace3 - 2*ybase(33)*10^6);
plot3(x, g1m13_medconc*ones(1, length(x))*10^6, syn_trace1, 'LineStyle','-', 'LineWidth',4,'color','w');
plot3(x, g1m3_medconc*ones(1, length(x))*10^6, syn_trace2, 'LineStyle','-', 'LineWidth',4,'color','k');
plot3(x, g1m1_medconc*ones(1, length(x))*10^6, syn_trace3, 'LineStyle','-', 'LineWidth',4,'Color',[128 128 128]/256);
% Calculate and plot the synergy 'slice' at each affinity level
syn_traceg18 = (trace_g18 - ybase(33)*10^6) - (add_traceg18 - 2*ybase(33)*10^6);
syn_traceg7 = (trace_g7 - ybase(33)*10^6) - (add_traceg7 - 2*ybase(33)*10^6);
plot3(glyc_kon{"FcgRIIIA-158V", "G18"}*ones(1, length(x))*10^-6, y, syn_traceg18, 'LineStyle','-', 'LineWidth',4)
plot3(glyc_kon{"FcgRIIIA-158V", "G2"}*ones(1, length(x))*10^-6, y, syn_traceg7, 'LineStyle','-', 'LineWidth',4,'color',[0 114 188]/256)
legend(["", "G1M1,3", "G1M3", "G1M1", "High affinity glycosylation", "Baseline Fc\gammaRIIIa-V^{158}"]);
title("Difference between A and B")
view(-20,35)
hold off;
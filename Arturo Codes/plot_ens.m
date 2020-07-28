function plot_ens(ensmat_in, enscells_in, ncellsperens, legends)
% Plots the activation sequence and core cells of ensembles.
% INPUT:
%   ensmat_in: <logical> nens-by-T-by-M matrix with the activity of each ensemble.
%   enscells_in: <logical> N-by-nens-by-M matrix with the the core cells of each ensemble.
%   ncellsperens: <integer> Nens-by-1 vector with the number of core cells per ensemble.
%     If it's a scalar then all ensembles have the same number of core cells. If it's equal to 0, 
%     the number is generated randomly. 
%   legends: [optional] <string> 1-by-M cell with the legends to add in the plot.

% Define variables
N = size(enscells_in, 1); %Number of neurons in population
T = size(ensmat_in, 2); %Number of total frames
nens = size(ensmat_in, 1); %Number of ensembles
M = size(ensmat_in, 3); %Number of different activation sequences/core cells detected (1: real, 2: carrillo, 3: herzog)

info_str = "| {\\color[rgb]{.5 .5 0}N = %d} | {\\color[rgb]{.5 0 .5}Nens = %d} | {\\color[rgb]{0 .5 .5}Ncorecells = %s} |";
plot_info = sprintf(info_str,N,nens,mat2str(ncellsperens));
plot_title1 = "Ensemble Activacion Sequence";
plot_title2 = "Core Neurons per Ensemble";

colors = {'b','r','g','y','m'};

%% Plot Ensemble Activation Sequence 
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for m = 1:M
    delta = ((-1)^(m-1))*0.2*round((m-1)/2);
    [y_axis, x_axis] = find(ensmat_in(:,:,m) == 1);
    scatter(x_axis, y_axis+delta, colors{m}, 'square', 'filled');
end

% Configuration of plot
title({plot_title1,plot_info});
xlabel('# Frame'); ylabel('# ID Ensemble');
xlim([0,T]); ylim([0,nens+0.5]);
if exist('legends','var')
    legend(legends);
end

% Plot Core Neurons
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for i = 1:nens
    for m = 1:M
        x_axis = find(enscells_in(:,i,m) == 1);
        y_axis = i*ones(length(x_axis),1);
        delta = ((-1)^(m-1))*0.2*round((m-1)/2);
        scatter(x_axis, y_axis+delta, colors{m}, 'square', 'filled');
    end
end

% Configuration of plot
title({plot_title2,plot_info});
xlabel('# Neuron'); ylabel('# ID Ensemble');
xlim([0,N]); ylim([0,nens+0.5]);
if exist('legends','var')
    legend(legends);
end

end
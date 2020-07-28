clear;clc;close all;
%% Parameters for raster
T = 1000; %Number of frames
fr = 0.2; %Firing rate of neurons
nens = 12; %Number of ensembles
time_inactivity = 0.2; %Percentage of time without ensemble
ncellsperens = 35; %Number of neurons per ensemble
N_min = 300; %Minimum value of neuron population
N_max = 800; %Maximum value of neuron population
N_size = 2; %Number of different values of neuron population
N = linspace(N_min,N_max,N_size)'; %Vector with number of neurons population (variable)
R = 2; %Number of repetitions

%% Pre-allocated memory
n_metrics = 5; %Number of metrics to calculate
M = 2; %Number of algorithms to compare
performance_fpr_tpr = zeros(N_size,n_metrics,R,M); %Matrix with performance metrics using TPR and FPR
performance_jaccard = zeros(N_size,n_metrics,R,M); %Matrix with performance metrics using Jaccard distance
nens_alg = zeros(M,1); %Vector with number of ensembles detected by algorithms
times_comp = zeros(M,1); %Vector with computing time of algorithms
ensmat_in = false(nens,T,M+1); %Matrix with activation sequence of each ensemble

t_init = tic;
for i = 1:N_size
    t_0 = tic;
    for r = 1:R
        %% Local variables inicialization
        enscells_in = false(N(i),nens,3);
        %% Generate raster
        [ensmat_in(:,:,1),enscells_in(:,:,1),raster,frates] = generate_data(N(i),T,nens,ncellsperens,fr,time_inactivity);
        %% Ensembles detection and sorting
        [ensmat_in_alg, enscells_in_alg, nens_alg(1), times_comp(1)] = get_carrillo_ens(raster);
        [ensmat_in(:,:,2),enscells_in(:,:,2)] = sort_by_equivalence(ensmat_in(:,:,1),enscells_in(:,:,1),ensmat_in_alg, enscells_in_alg,'cosine');
        [ensmat_in_alg, enscells_in_alg, nens_alg(2), times_comp(2)] = get_herzog_ens(raster);
        [ensmat_in(:,:,3),enscells_in(:,:,3)] = sort_by_equivalence(ensmat_in(:,:,1),enscells_in(:,:,1),ensmat_in_alg, enscells_in_alg,'cosine');
        %% Get performance
        % Using FPR and TPR metric
        performance_alg = get_performance(ensmat_in,enscells_in,'fpr_tpr');
        performance_fpr_tpr(i,:,r,1) = [nens_alg(1),times_comp(1),performance_alg(1,:)];
        performance_fpr_tpr(i,:,r,2) = [nens_alg(2),times_comp(2),performance_alg(2,:)];
        % Using Jaccard metric
        performance_alg = get_performance(ensmat_in,enscells_in,'jaccard');
        performance_jaccard(i,:,r,1) = [nens_alg(1),times_comp(1),performance_alg(1,:)];
        performance_jaccard(i,:,r,2) = [nens_alg(2),times_comp(2),performance_alg(2,:)];
    end
    t_i = toc(t_0);
    fprintf('Iteration %d (%d seconds)\n',i,t_i);
end
t_elapsed = toc(t_init);
fprintf('Total time elapsed: %d seconds\n',t_elapsed);

%% Plot performance curves
legends = {'Carrillo', 'Herzog'};
x_label = 'N (neurons)';
F_str = sprintf('T = %d',T);
% Plot number of ensembles detected and computing time
plot_performance(N,performance_fpr_tpr(:,1,:,:),x_label,'Detected Ensembles','Ensemble Detection',legends,F_str,nens);
plot_performance(N,performance_fpr_tpr(:,2,:,:),x_label,'Computing Time [seconds]','Ensemble Detection',legends,F_str);
% Plot AUC of global activation, individual activation and core-cells
plot_performance(N,performance_fpr_tpr(:,3,:,:),x_label,'AUC','Global Activation Times Detection',legends,F_str);
plot_performance(N,performance_fpr_tpr(:,4,:,:),x_label,'AUC','Individual Activation Times Detection',legends,F_str);
plot_performance(N,performance_fpr_tpr(:,5,:,:),x_label,'AUC','Core-Cells Detection',legends,F_str);
% Plot Jaccard coefficient of global activation, individial activation and core-cells
plot_performance(N,performance_jaccard(:,3,:,:),x_label,'Mean of Jaccard Similarity','Global Activation Times Detection',legends,F_str);
plot_performance(N,performance_jaccard(:,4,:,:),x_label,'Mean of Jaccard Similarity','Individual Activation Times Detection',legends,F_str);
plot_performance(N,performance_jaccard(:,5,:,:),x_label,'Mean of Jaccard Similarity','Core-Cells Detection',legends,F_str);
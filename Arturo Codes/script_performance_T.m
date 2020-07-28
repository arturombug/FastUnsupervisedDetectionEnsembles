clear;clc;close all;
%% Parameters for raster
N = 300; %Number of neurons
fr = 0.2; %Firing rate of neurons
nens = 12; %Number of ensembles
time_inactivity = 0.2; %Percentage of time without ensemble
ncellsperens = 35; %Number of neurons per ensemble
T_min = 1000; %Minimum value of frames
T_max = 1500; %Maximum value of frames
T_size = 5; %Number of different values of frames
T = linspace(T_min,T_max,T_size); %Vector with number of frames (variable)
R = 3; %Number of repetitions

%% Pre-allocated memory
n_metrics = 5; %Number of metrics to calculate
M = 2; %Number of algorithms to compare
performance_fpr_tpr = zeros(T_size,n_metrics,R,M); %Matrix with performance metrics using TPR and FPR
performance_jaccard = zeros(T_size,n_metrics,R,M); %Matrix with performance metrics using Jaccard distance
nens_alg = zeros(M,1); %Vector with number of ensembles detected by algorithms
times_comp = zeros(M,1); %Vector with computing time of algorithms
enscells_in = false(N,nens,M+1); %Matrix with core-cells of each ensemble

t_init = tic;
for i = 1:T_size
    t_0 = tic;
    for r = 1:R
        %% Local variables inicialization
        ensmat_in = false(nens,T(i),3);
        %% Generate raster
        [ensmat_in(:,:,1),enscells_in(:,:,1),raster,frates] = generate_data(N,T(i),nens,ncellsperens,fr,time_inactivity);
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

%% Plot performances curves
legends = {'Carrillo', 'Herzog'};
x_label = 'T (bins)';
F_str = sprintf('N = %d',N);
% Plot number of ensembles detected and computing time
plot_performance(T,performance_fpr_tpr(:,1,:,:),x_label,'Detected Ensembles','Ensemble Detection',legends,F_str,nens);
plot_performance(T,performance_fpr_tpr(:,2,:,:),x_label,'Computing Time [seconds]','Ensemble Detection',legends,F_str);
% Plot AUC of global activation, individual activation and core-cells
plot_performance(T,performance_fpr_tpr(:,3,:,:),x_label,'AUC','Global Activation Times Detection',legends,F_str);
plot_performance(T,performance_fpr_tpr(:,4,:,:),x_label,'AUC','Individual Activation Times Detection',legends,F_str);
plot_performance(T,performance_fpr_tpr(:,5,:,:),x_label,'AUC','Core-Cells Detection',legends,F_str);
% Plot Jaccard coefficient of global activation, individial activation and core-cells
plot_performance(T,performance_jaccard(:,3,:,:),x_label,'Mean of Jaccard Similarity','Global Activation Times Detection',legends,F_str);
plot_performance(T,performance_jaccard(:,4,:,:),x_label,'Mean of Jaccard Similarity','Individual Activation Times Detection',legends,F_str);
plot_performance(T,performance_jaccard(:,5,:,:),x_label,'Mean of Jaccard Similarity','Core-Cells Detection',legends,F_str);
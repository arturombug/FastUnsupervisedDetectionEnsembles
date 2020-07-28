clear;clc;close all;
%% Parameters for raster
N = 300; %number of neurons in the population
fr = 0.2; %firing rate of neurons
nens = 12; %number of ensembles
time_inactivity = 0.2;
ncellsperens = 35;
T = linspace(1000,1500,5);
T_size = length(T);

%% Pre-allocated memory
performance_fpr_tpr = zeros(T_size,5,2);
performance_jaccard = zeros(T_size,5,2);
nens_alg = zeros(2,1);
times_comp = zeros(2,1);
enscells_in = false(N,nens,3);
t_init = tic;
for i = 1:T_size
    t_0 = tic;
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
    performance_alg = get_performance(ensmat_in,enscells_in,'fpr_tpr');
    performance_fpr_tpr(i,:,1) = [nens_alg(1),times_comp(1),performance_alg(1,:)];
    performance_fpr_tpr(i,:,2) = [nens_alg(2),times_comp(2),performance_alg(2,:)];
    performance_alg = get_performance(ensmat_in,enscells_in,'jaccard');
    performance_jaccard(i,:,1) = [nens_alg(1),times_comp(1),performance_alg(1,:)];
    performance_jaccard(i,:,2) = [nens_alg(2),times_comp(2),performance_alg(2,:)];
    t_i = toc(t_0);
    fprintf('Iteration %d (%d seconds)\n',i,t_i);
end
t_elapsed = toc(t_init);
fprintf('Total time elapsed: %d seconds',t_elapsed);

%% Plot performances curves
legends_0 = {'Carrillo', 'Herzog', 'Real'};
legends = {'Carrillo', 'Herzog'};
x_label = 'T (bins)';
plot_performance(T,[performance_fpr_tpr(:,1,1),performance_fpr_tpr(:,1,2)]',x_label,'Detected Ensembles','Ensemble Detection',legends,nens);
plot_performance(T,[performance_fpr_tpr(:,2,1),performance_fpr_tpr(:,2,2)]',x_label,'Computing Time [seconds]','Ensemble Detection',legends);

plot_performance(T,[performance_fpr_tpr(:,3,1),performance_fpr_tpr(:,3,2)]',x_label,'AUC','Global Activation Times Detection',legends);
plot_performance(T,[performance_fpr_tpr(:,4,1),performance_fpr_tpr(:,4,2)]',x_label,'AUC','Individual Activation Times Detection',legends);
plot_performance(T,[performance_fpr_tpr(:,5,1),performance_fpr_tpr(:,5,2)]',x_label,'AUC','Core-Cells Detection',legends);

plot_performance(T,[performance_jaccard(:,3,1),performance_jaccard(:,3,2)]',x_label,'Mean of Jaccard Distance','Global Activation Times Detection',legends);
plot_performance(T,[performance_jaccard(:,4,1),performance_jaccard(:,4,2)]',x_label,'Mean of Jaccard Distance','Individual Activation Times Detection',legends);
plot_performance(T,[performance_jaccard(:,5,1),performance_jaccard(:,5,2)]',x_label,'Mean of Jaccard Distance','Core-Cells Detection',legends);
function [perfdata, enscells_in, enscells_out] = SyntEnsTest_DensityPeaks_Vdc(N,T,nens,fr,ncellsperens,ntimesperens,dc,nsur,prct)
%
% [perfdata] = SyntEnsTest_DensityPeaks(N,T,nens,fr,ncellsperens,ntimesperens,dc,thr)
% Generates synthetic ensembles and then detects the ensembles using search
% of density peaks clustering method.
%
% INPUTS
%
% 'ncells' is the number of cells of the raster
% 'fr' is either a ncells x 1 vector with the neuron firing rates or a
% scalar representing the mean firing rate of the population.
% 'T' is the duration in bins of the raster
% 'nens' is the number of ensembles to plug in the raster
% 'ncellsperens' is a nens x 1 vector with the number of cells per
% ensemble.
% 'ntimesperens' is a nens x 1 vector with the proportion [0,1] of the recording
% that each ensemble is active.
% 'dc' average percentage of neighbours, ranging from [0, 1], must be a
% vector
% 'nsur' is the number of surrogates to compute core-cells
% 'prct' is a vector or scalar with the percentiles used as threshold on
% surrogate data
%
%
% OUTPUT
%
% 'perfdata' structure with TPR and FPR for ensemble activation times and
% core cells.


% 1.- Generating ensembles
[ensmat_in,enscells_in,raster] =  MakeEnsembles_fix_rate(N,fr,T,nens,ncellsperens,ntimesperens); % genrates synthetic ensembles

% variables to store data
ndc = length(dc);
nprct = length(prct);
tpr_times = zeros(nens,length(ndc));
fpr_times = tpr_times;
out_ens = zeros(ndc,1);

tpr_cells = zeros(nens,nprct);
fpr_cells = tpr_cells;

% 2.- Detecting core-cells using correlation and surrogate data
[enscells_out] = core_cells_by_corr(raster,ensmat_in,nsur,prct);

for n=1:nens % loop over ensembles
    t_cells = find(enscells_in(:,n)); % indices of input ensemble core cells
    for p=1:nprct % loop over thresholds
        e_cells = find(enscells_out(:,n,p));% indices of output ensemble core cells
        [tpr_cells(n,p),fpr_cells(n,p)] = tpr_fpr_pks(t_cells,e_cells); % computing tpr and fpr
    end
end


% 3.- Clustering
[~,score] = pca(raster','numcomponent',6); % pca with 6 components
distM = pdist2(score,score); % distance matrix based on pca

for i=1:length(dc)    
    [Nens,ensId]=find_ens_by_density(distM,dc(i),-1,0); % clustering
    out_ens(i)=Nens;
    
    % Computes correlation between input and output ensemble sequences
    % Matches output ensembles to the input ensembles based on maximal
    % correlation between input and output ensembles
    ensmat_out = bsxfun(@eq,ensId,(1:Nens)');
    % checking that input and output have the same shape
    if size(ensmat_out,2)~=size(ensmat_in,2);
        ensmat_out = ensmat_out';
    end
    C = 1-pdist2(single(ensmat_in),single(ensmat_out),'correlation');
    [~,ens_id] = max(C,[],2);
    
    % If any, here we look for the output ensembles that were not plugged into the raster
    extraens =[];
    if Nens>nens
        ee = 1:Nens;
        extraens = find(~ismember(ee,ens_id));
    end
     
    % 4.- Performance: Ture Positive and False Positive Rates
    % for ensemble activation times         
    for n=1:nens
        t_times = find(ensmat_in(n,:)); % indices of input ensemble activations times
        e_times = find(ensmat_out(ens_id(n),:)); % indices of output ensemble activations times
        [tpr_times(n,i),fpr_times(n,i)] = tpr_fpr_pks(t_times,e_times); % computing tpr and fpr
    end

end

% storing data on perfdata
perfdata.tpr_times = tpr_times;
perfdata.fpr_times = fpr_times;
perfdata.tpr_cells = tpr_cells;
perfdata.fpr_cells = fpr_cells;
perfdata.prct = prct;
perfdata.dc = dc;
perfdata.input_ens = nens;
perfdata.output_ens = out_ens;
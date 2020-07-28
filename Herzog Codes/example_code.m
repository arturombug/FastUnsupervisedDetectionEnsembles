% EXAMPLE CODE FOR USING neuronRaster2EnsRaster
% Rubén Herzog A 2017.

% loading sample data
filename = 'natural_im_data';
load([filename,'.mat'])


% setting parameters
pars.sr = 1; % Hz
pars.bin = 0.005; % seconds
pars.dc = 0.01; % cut-off for distances
pars.ishalo = 0; % not using halo refinement of the clusters
pars.npcs = 6;
pars.maxmem = 4*10^9; % maximal 4Gb of RAM 
pars.minspk = 3; % minimum 3 spikes per pattern
pars.nsur = 10; % surrogates for core-cells, should be 1000 or more.
pars.prct = 99;% percentile on the surrogate core-cell distribution
pars.ccthr = 0.1; % minimal correlation between template and pattern
pars.sampfac = 50; % sampling factor; controls how many subsamples ara drawn.
pars.thrmet = 'fit'; % method to automatically detect centroids

% running analysis
[ens_seq,maxcor,templates,core_cells,clust_out,pars] = neuronRaster2EnsRaster(ts,filename,pars);
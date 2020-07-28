function [ensmat_in,enscells_in,rasters_fr,frates,st_pat,ens_ras] = Make_ST_Ensembles_set_of_frates(ncells,fr,T,nens,ncellsperens,ntimesperens,ens_r)
%
% [ensmat_in,enscells_in,raster,frates] = Make_ST_Ensembles_fix_rate(ncells,fr,T,nens,ncellsperens,ntimesperens,ens_r)
% Generates a raster with random ensemble activity. Each ensemble has range
% 'ens_r', i.e., the number of bins the ensemble is extended in time. The
% firing rate distribution is drawn from a rectified gaussian distribution
% of zero mean and std fr.
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
% 'ens_r' is a nens x 1 vector, where 'ens_r(i)' is the range of the
% ensemble i. The range is the number of bins the pattern uses.
%
% OUTPUTS
% 
% 'ensmat_in' is nens x T  binary matrix with ground truth ensembles
% activation
% 'enscells_in' is a ncells x nens binary matrix with (i,j)=1 if neuron i
% participates on ensemble j.
% 'raster' is the generated raster with the ensembles
if sum(ntimesperens.*ens_r)>1
    error('sum(ntimesperens*ens_r) must be lower than 1')
    
elseif any(ncellsperens>ncells)
    error('ncellsperens be lower than ncells')
    
        
end

% 1.- Generating the firing rates vector
tt = 1:T; % timebase of the raster
nfrates = length(fr);
frates = zeros(ncells,nfrates);
if nfrates== ncells % fr is a vector with frs
    aux_frates = reshape(fr,[ncells 1]); % just to be sure that has the right shape
    
else % fr is scalar or set of frates
    for f=1:nfrates
        aux_frates = 1;
        while length(aux_frates)<ncells % looping while the number of positive entries of frates is lower than ncells
            aux_frates = fr(f)*randn(ncells*10,1); % normally distributed firing rates with mean avefr
            aux_frates(aux_frates<0.001 | aux_frates>0.5 ) =[]; % removing negative firing rates and very high firing rates
        end
        frates(:,f) = aux_frates(randperm(length(aux_frates),ncells)); % extracting random values from the positive part of the distribution
    end
end

% 2.- Generating ensembles (cells participating), activation time and plugging ensembles on raster
st_pat = cell(nens,1);
raster = false(ncells,T);
enscells_in = false(ncells,nens); % binary matrix with (i,j)=1 if neuron i participates on ensemble j
ensmat_in = false(nens,T); % binary matrix with (i,j)=1 if ensemble i was active on the significant bin 
tt2 = tt;
tmax = T;
for i=1:nens
        % 2.1 enerating the ensemble cells        
        ens_patt = false(ncells,ens_r(i));     
        cells = [];
        for r=1:ens_r(i)
            cells = cat(2,cells,(randperm(ncells,floor(ncellsperens(i)./ens_r(i))))'); 
            ens_patt(cells(:,r),r) = true;
        end
        enscells_in(cells(:),i) = true;
        st_pat{i} = ens_patt;
        % 2.2 generating the ensemble activation times        
        %enspos = randperm(tmax-ens_r(i)-1,floor(ntimesperens(i)*T)); % generating random indices for the ensembles activation
        [enspos] = randi_with_int(floor(ntimesperens(i)*T),tmax-ens_r(i)-1,ens_r(i));       

        % 2.3 Plugging the ensemble on the raster     
        rem_tt=[];
        for r=1:ens_r(i)
            ensmat_in(i,tt2(enspos+r-1)) = true; % assigning to the activation binary matrix
            raster(:,tt2(enspos+r-1)) = repmat(ens_patt(:,r),[1 length(tt2(enspos+r-1))]); % plugging the ensemble
            rem_tt = cat(2,rem_tt,enspos+r-1);
        end
        tt2(rem_tt(:))=[]; % removing those indices from the timebase to avoid repetitions        
        tmax = length(tt2); % maximum number of indices after removal
    
end
ens_ras = raster;
rasters_fr = false(ncells,T,nfrates);
fr_ens = sum(raster,2);
% 3.- Checking excess or deficit of spikes
for f=1:nfrates    
    E_fr = floor(frates(:,f).*T);
    fr_dif = fr_ens - E_fr; % negative is excess and positive deficit
    rem_spk = find(fr_dif>0); % positive means spike removal
    add_spk = find(fr_dif<0); % negative means spike addition
    rasters_fr(:,:,f) = raster;
    
    % 3.2 Adding spikes
    for n=1:length(add_spk)
        no_spks = find(raster(add_spk(n),:)==0); % bins without spikes
        add_times = randperm(length(no_spks),abs(fr_dif(add_spk(n)))); % abs(fr_dif) random indices
        rasters_fr(add_spk(n),no_spks(add_times),f)=true;
    end
    
    % 3.3 removing spikes
    for n=1:length(rem_spk)
        yes_spks = find(raster(rem_spk(n),:)); % bins with spikes
        rem_times = randperm(length(yes_spks),abs(fr_dif(rem_spk(n)))); % abs(fr_dif) random indices
        rasters_fr(rem_spk(n),yes_spks(rem_times),f)=false;
    end
end



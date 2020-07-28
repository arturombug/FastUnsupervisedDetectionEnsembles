%% Detecting spatio-temporal patterns with clustering
% the idea is to generate a raster with patterns with range >=1. Then, the
% patterns will be detected clustering bins as if they were independent.
% Finally, based on the correlation between the ensembles, the
% spatio-temporal sequence if built-up from the correlated ensembles.
%
%% 1.- Raster variables
N=100;
T=5000;
fr=0.2;
nreps = 20;
ntimesperens = [0.066 0.1 0.2 0.2]; % probability of each ensemble
ncellsperens = [40 20 10 10]; % cells per ensemble
ens_r = [3 2 1 1];
dc = 0.008:0.002:0.06;
npcs = [3 4 5 6 7 8];
ndc = length(dc);
nnpcs = length(npcs);

%% running
for r=1:nreps
    [ensmat_in,enscells_in,raster,frates_fix,st_patt] = Make_ST_Ensembles_fix_rate(N,fr,T,nens,ncellsperens,ntimesperens,ens_r); % genrates synthetic ensembles
    selbins = sum(raster)>minspk;
    truelabels= sum(bsxfun(@times,ensmat_in(:,selbins)',1:nens),2);
    raster_sel = raster(:,selbins);
    ensmat_in = ensmat_in*1;
    for p=1:nnpcs
        [~,pcs] = pca(raster_sel','NumComponents',npcs(p)); % pca with npcs num components
        bincor = pdist2(pcs,pcs); % euclidean distance on principal component space
        parfor d=1:ndc;
            [~, rho] = paraSetv2(bincor, dc(d) );
            delta = delta_from_dist_mat(bincor, rho);
            % 2.4 clustering for each set of delta and rho
            [numClust,centInd,predbounds] = cluster_by_pow_fit(delta,rho,pb);
            % each point is assigned to its closest centroid
            if numClust==1
                labels = ones(length(delta),1);
            else
                dist2cent = bincor(centInd>0,:); % distance from centroid to any other point
                [~,labels] = min(dist2cent);
            end
            cents{d,p,r} = arrayfun(@(x) find(centInd==x),1:numClust);
            
            % ensemble raster
            ensmat_out = zeros(numClust,T);
            ensmat_out(:,selbins) = bsxfun(@eq,labels',(1:numClust))';
            
            % Detecting core cells
            [core_cells,~,~] = find_core_cells_by_correlation(raster,ensmat_out,nsur,prct);
            id_sel_core = sum(core_cells,1)>minspk;
            cents{d,p,r} = cents{d,p,r}(id_sel_core);
            
            % re-assigning cluster to fit the input clusters
            ensmat_out = ensmat_out(id_sel_core,:);
            C = 1-pdist2(ensmat_in,ensmat_out,'correlation');
            [~,ens_id] = max(C,[],2);
            % If any, here we look for the output ensembles that were not plugged into the raster
            if numClust>nens
                ee = 1:numClust;
                extraens = find(~ismember(ee,ens_id));
            else
                extraens =[];
            end
            % validation
            labels_sel = sum(bsxfun(@times,ensmat_out(:,selbins),(1:size(ensmat_out,1))'));
            
            meanErr(d,p,r) = 1 - mean(sum(ensmat_in == ensmat_out(ens_id,:),2)./length(ensmat_in));
            clus_validation(d,:,p,r) =  cluster_validation_indexes(pcs,bincor,labels_sel,cents{d,p,r},length(cents{d,p,r}),'euclidean');
            labels_all{d,p,r} = labels_sel;
        end
    end
end
function [numClust, centInd,predbounds] = cluster_by_pow_fit(delta,rho,pb)
% fits a power law to the rho vs delta plot and uses 'pb' as a threshold
% for detecting cluster centroids.
NE = length(rho);
centInd = zeros(1, NE);
% fitting a line to log(delta) vs log(rho) and using 99.9
% upper confidence interval of fit as threshold

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'bisquare';
mindelta = 10^-4; % min delta to be considered, improves the fit
nzind = find(delta>mindelta & rho>0); % uses only points with delta>mindelta and rho>0
nzdelta = double(delta(nzind));% just to make sure they are double
nzrho = double(rho(nzind));

% fitting in log space
[fitresult] = fit(log(nzrho)',log(nzdelta)', ft, opts );
predbounds = predint(fitresult,log(nzrho),pb/100,'observation','on'); % prediction bound computation

% selecting as centroids points above the confidence bound,
auxid = (nzdelta>exp(predbounds(:,2))');
predbounds = cat(2,nzrho',exp(predbounds(:,2)));
centid = nzind(auxid); % indices on the original basis
numClust = length(centid);
centInd(centid) = 1:numClust;
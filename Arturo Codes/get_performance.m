function performance = get_performance(ensmat_in,enscells_in,method,plot_titles,legends)
% Gets the performance metrics.
% INPUT:
%   ensmat_in: <logical> nens-by-T-by-(M+1) matrix with the activity of each ensemble.
%     M is the number of algorithms (there is one more that corresponds to the real values).
%   enscells_in: <logical> N-by-nens-by-(M+1) matrix with the core cells of each ensemble.
%     M is the number of algorithms (there is one more that corresponds to the real values).
%   method: <string> with the name of method used to calculate the performance. Available
%     options are 'fpr_tpr' and 'jaccard'.
%   plot_title: [optional] <string> with the title of plot (only for TPR-FPR option).
%   legends: [optional] <string> 1-by-M cell with the legends to add in the plot (only for TPR-FPR option).

% Get the number of algorithms to evaluate and pre-allocate memory
M = size(ensmat_in,3)-1;
performance = zeros(M,3);
nens = size(enscells_in,2);

% Get performance metrics
ensmat_global = logical(sum(ensmat_in,1));
switch method
    % If TPR-FPR method is used
    case 'fpr_tpr'
        % Pre-allocate memory
        fpr_global = zeros(M,1);
        tpr_global = zeros(M,1);
        fpr_indiv = zeros(M,nens);
        tpr_indiv = zeros(M,nens);
        fpr_core = zeros(M,nens);
        tpr_core = zeros(M,nens);
        for i = 1:M
            isens = any(ensmat_in(:,:,i+1),2); %Vector indicating which ensembles were detected
            % Calculate TPR-FPR for global activation sequence
            [fpr_global(i,:), tpr_global(i,:)] = fpr_tpr(ensmat_global(:,:,1)', ensmat_global(:,:,i+1)');
            % Calculate TPR-FPR for individual activation sequence
            [fpr_indiv(i,isens'), tpr_indiv(i,isens')] = fpr_tpr(ensmat_in(isens,:,1)', ensmat_in(isens,:,i+1)');
            % Calculate TPR-FPR for core-cells
            [fpr_core(i,isens'), tpr_core(i,isens')] = fpr_tpr(enscells_in(:,isens',1), enscells_in(:,isens',i+1));
        end
        % If plot_titles variable exists, then plot ROC curves
        if exist('plot_titles','var')
            AUC_sg = plot_roc(fpr_global, tpr_global, plot_titles{1}, legends);
            AUC_si = plot_roc(sort(fpr_indiv,2), sort(tpr_indiv,2), plot_titles{2}, legends);
            AUC_cc = plot_roc(sort(fpr_core,2), sort(tpr_core,2), plot_titles{3}, legends);
        else
            AUC_sg = plot_roc(fpr_global, tpr_global);
            AUC_si = plot_roc(sort(fpr_indiv,2), sort(tpr_indiv,2));
            AUC_cc = plot_roc(sort(fpr_core,2), sort(tpr_core,2));
        end
        performance = [AUC_sg, AUC_si, AUC_cc];
    % If Jaccard distance is used
    case 'jaccard'
        for i = 1:M
            isens = any(ensmat_in(:,:,i+1),2); %Vector indicating which ensembles were detected
            % Calculate Jaccard similarity for global activation sequence
            jaccard_sg = 1-pdist2(double(ensmat_global(:,:,1)),double(ensmat_global(:,:,i+1)),'jaccard');
            % Calculate Jaccard similarity for individual activation sequence
            jaccard_si = mean(diag(1-pdist2(double(ensmat_in(isens,:,1)),double(ensmat_in(isens,:,i+1)),'jaccard')));
            % Calculate Jaccard similarity for core-cells
            jaccard_cc = mean(diag(1-pdist2(double(enscells_in(:,isens',1)'),double(enscells_in(:,isens',i+1)'),'jaccard')));
            performance(i,:) = [jaccard_sg, jaccard_si, jaccard_cc];
        end
end

end
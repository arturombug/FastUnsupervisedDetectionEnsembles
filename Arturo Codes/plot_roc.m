function AUC = plot_roc(fpr, tpr, plot_title, legends)
% Plot the ROC curve and return the Area Under Curve (AUC).
% INPUT:
%   fpr: <float> N-by-M vector with FPR (False Positive Rate). N is the number of data sets and M
%     is the length of data sets.
%   tpr: <float> N-by-M vector with TPR (True Positive Rate). N is the number of data sets and M
%     is the length of data sets.
%   plot_title: [optional] <string> with the title of plot.
%   legends: [optional] <string> 1-by-M cell with the legends to add in the plot.
% OUTPUT:
%   AUC: <float> N-by-1 vector with area under ROC curve of each data set.

% Get length of different data
N = size(fpr,1);
% Pre-allocate memory
AUC = zeros(N,1);

% If you want to plot ROC curve
if exist('legends','var')
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    % Plot ROC curve for each data set in the same figure
    for i = 1:N
        x_axis = [0, fpr(i,:), 1];
        y_axis = [0, tpr(i,:), 1];
        plot(x_axis, y_axis, 'LineWidth', 1.5);
        % Get AUC
        cumulative_area = cumtrapz(x_axis, y_axis);
        AUC(i) = cumulative_area(end);
    end
    % Configuration of plot
    title(plot_title);
    xlabel('FPR');
    ylabel('TPR');
    xlim([0,1]);
    ylim([0,1]);
    legend(legends);    
% If you don't want to plot ROC curve
else 
   for i = 1:N
       % Get AUC
       x_axis = [0, fpr(i,:), 1];
       y_axis = [0, tpr(i,:), 1];
       cumulative_area = cumtrapz(x_axis, y_axis);
       AUC(i) = cumulative_area(end);
   end
end

end
function plot_performance(X,Y,variable_str,performance_str,plot_title,legends,F_str,real_value)
% Plots curves of performance.
% INPUT:
%   X: <integer> X_size-by-1 vector with different values of X variable.
%   Y: <integer> X_size-by-R-by-M matrix with performance metric for each value of X and
%     each repetition. M is the number of algorithms to evaluate.
%   variable_str: <string> with the name of the variable parameter.
%   performance_str: <string> with the name of the metric to plot.
%   plot_title: <string> with the title of plot.
%   legends: <string> 1-by-M cell with the legends to add in the plot.
%   F_str: <string> with the fixed parameter and its value, e.g. N = 100.
%   real_value: [optional] <float> with the real value of the metric.

% Initialized variables
Y = squeeze(Y);
if size(Y,3) > 1
    R = size(Y,2);
    M = size(Y,3);
else
    R = 1;
    M = size(Y,2);
end
colors = {'b','r','g','y','m'};
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;

% Calculate media and standard deviation of data
if R > 1
    Y_mean = mean(Y,2);
    Y_std = std(Y,0,2);
end
% Plot mean and confidence interval of metric value versus variable parameter
for m = 1:M
    if R > 1
        P(m) = plot(X,Y_mean(:,m),'Marker','.','MarkerSize',20,'Color',colors{m});
        A = area(X,[Y_mean(:,m)-4*Y_std(:,m),8*Y_std(:,m)],'LineStyle','-.','LineWidth',1.5,'FaceColor',colors{m},'EdgeColor',colors{m});
        A(1).FaceAlpha = 0; A(2).FaceAlpha = 0.15; 
        A(1).EdgeAlpha = 0; A(2).EdgeAlpha = 0.2; 
    else
        P(m) = plot(X,Y(:,m),'Marker','.','MarkerSize',20,'Color',colors{m});
    end
end

% Configuration of plot
title(sprintf('%s {\\color[rgb]{.5 .5 0}(%s)}',plot_title,F_str),'FontSize',13);
xlabel(variable_str);ylabel(performance_str);
if R > 1
    ylim([1.05*min(0,min(Y_mean(:)-4*(Y_std(:)))),1.05*max(Y_mean(:)+4*Y_std(:))]);
else
    ylim([0,1.02*max(Y(:))]);
end
xlim([0.98*min(X),1.02*max(X)]);

% Plot the real value if it exists
if exist('real_value','var')
    plot(X,real_value*ones(length(X),1),'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);
end

% Add legends to plot
legend(P,legends); 

end
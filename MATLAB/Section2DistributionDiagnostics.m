%% Severity Analysis: Graphical Diagnostics
% Data: 25 non-zero casualty events (1916-2022)
% Diagnostics: Mean Excess, Max-to-Sum, Hill Estimator (user script), QQ Plot vs. Normal
% References: Embrechts et al. (1997), Cirillo & Taleb (2016)
% Mean Excess per Cirillo & Taleb (2016, p. 30): e(u) = E[max(Y - u, 0)] / P(Y > u)
close all
clear

% Load and prepare data
load('ANdata.mat')
nonzero_casualties = dvec(dvec(:,2) > 0, 2); % 25 non-zero casualties
xraw = nonzero_casualties; % Use raw data since theta=0 is fitted

% Mean Excess Plot (Cirillo & Taleb (2016) definition)
sorted_data = sort(xraw);
u_range = sorted_data(1:end-1); % Possible thresholds u < max
mex = zeros(size(u_range));
for i = 1:length(u_range)
    u = u_range(i);
    exceedances = xraw(xraw > u);
    mex(i) = mean(exceedances - u); % Corrected to match E[Y - u | Y > u], no extra division
end
figure;
t=tiledlayout(3,1)
nexttile(3)
plot(u_range, mex, 'b-o');
xlabel('Threshold u (Casualties)');
ylabel('Mean Excess e(u)');
title('Mean Excess Plot for Non-Zero Casualty Data');
%saveas(gcf, 'MeanExcessPlot.png', 'png');

% Max-to-Sum Plot
p_values = 1:4; % Moments 1 to 4
ms = zeros(length(xraw), length(p_values));
for p = 1:length(p_values)
    x = xraw .^ p_values(p);
    ms(:,p) = cummax(x) ./ cumsum(x);
end
%figure;

nexttile(2)
plot(ms);
xlabel('Number of Observations (n)');
ylabel('MS Ratio');
legend('Moment 1', 'Moment 2', 'Moment 3', 'Moment 4', 'Location', 'best');
legend('boxoff')
title('Max-to-Sum Plot for Non-Zero Casualty Data');
%saveas(gcf, 'MaxToSumPlot.png', 'png');

% Hill Estimator (using your script)
Hillestimator % call the script

% QQ Plot vs. Normal Distribution
pd_normal = fitdist(xraw, 'Normal');
theoretical_quantiles = norminv(((1:length(xraw)) - 0.5) / length(xraw), pd_normal.mu, pd_normal.sigma);
%figure;
figure
plot(theoretical_quantiles, sorted_data, 'bo'); hold on;
plot(theoretical_quantiles, theoretical_quantiles, 'r--');
xlabel('Theoretical Normal Quantiles');
ylabel('Empirical Quantiles');
title('QQ Plot: Non-Zero Casualties vs. Normal Distribution');
saveas(gcf, 'QQPlotNormal.png', 'png');
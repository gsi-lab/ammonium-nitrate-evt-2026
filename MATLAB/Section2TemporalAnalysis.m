%% Ammonium Nitrate Accident Analysis - Temporal Section
clear
close all
rng(1976); % Set random seed for reproducibility

%% Load and Prepare Data
T = readtable('MajorAccidents.xlsx', 'Sheet', 'Nitrate_2', 'Range', 'D1:H40'); % 39 rows, excluding header
T.Date = datetime(T.Date, 'InputFormat', 'MMMM d, uuuu'); % Parse full dates

% Compute decimal years (fractional)
start_dates = datetime(T.Year, 1, 1); % Year start
frac = (datenum(T.Date) - datenum(start_dates)) ./ (365 + (eomdate(T.Year, 2) == 29)); % Adjust for leap years
T.DecimalYear = T.Year + frac;

dvec = [T.DecimalYear, T.Deaths]; % Decimal year vs. casualties
[~, i2] = sort(dvec(:,1));
dvec = dvec(i2,:); % Sorted by decimal year

%% Inter-arrival Times for Non-Zero Casualty Events
nonzero_idx = find(dvec(:,2) > 0); % 25 non-zero casualty events
int_all = diff(dvec(nonzero_idx, 1)); % 23 inter-arrivals
% Note: First inter-arrival omitted by diff; minor effect with n=24

%% Fit Exponential Distribution and Estimate Parameters
pd = fitdist(int_all, 'Exponential'); % MLE fit
estmu = pd.mu; % Mean inter-arrival time
lambda = 1 / estmu; % Rate parameter

% Asymptotic 95% CI
ci = paramci(pd); % From covariance matrix


%% Dispersion Test for Homogeneity
decades = floor(dvec(nonzero_idx, 1) / 10) * 10; % Bin into decades
counts = histcounts(decades, [unique(decades); max(decades)+10]); % Counts per decade
mean_count = mean(counts);
var_count = var(counts);
dispersion_index = var_count / mean_count;
chi2_stat = (length(counts) - 1) * dispersion_index;
p_value = 1 - chi2cdf(chi2_stat, length(counts) - 1);
disp(['Dispersion Index: ', num2str(dispersion_index)]);
disp(['Chi-squared stat: ', num2str(chi2_stat), ', p-value: ', num2str(p_value)]);

%% Diagnostic Plots
figure;
tiledlayout(3,1);
ax = nexttile;
qqplot(ax, int_all, 'Exponential');
ylabel(ax, 'Quantiles for Inter-arrival Times');
nexttile;
t = linspace(0, max(int_all), 50); % Start at 0 for full range
histogram(int_all, 'BinMethod', 'auto', 'Normalization', 'pdf');
hold on;
plot(t, exppdf(t, estmu), '-');
xlabel('Interarrival Times (t, yr)');
ylabel('Probability Density, f(t)');
legend('Data', 'Exponential Fit', 'Location', 'best');
nexttile;
autocorr(int_all);
xlabel('Lag');
ylabel('ACF');
saveas(gcf, 'InterarrivalDiagnostics', 'png');

%% SAVE DATA
save ANdata dvec
%% supplementary material
% Bootstrap CI (10,000 iterations)
B = 10000;
bootlambdas = zeros(B, 1);
for i = 1:B
    bootdata = int_all(randi(length(int_all), length(int_all), 1));
    bootpd = fitdist(bootdata, 'Exponential');
    bootmu(i)=bootpd.mu;
    bootlambdas(i) = 1 / bootpd.mu;
end
ci_boot = prctile(bootlambdas, [2.5, 97.5]); % 95% bootstrap CI
ci_bootmu = prctile(bootmu, [2.5, 97.5]); % 95% bootstrap CI

% Histogram of Bootstrap Estimates for λ
figure;
histogram(bootlambdas, 50, 'EdgeColor', 'k');
xlabel('Bootstrap Estimates of λ');
ylabel('Frequency');
title('Histogram of Bootstrap Estimates for λ');
saveas(gcf, 'BootstrapLambdaHistogram', 'png');

%  POT Rate Parameter Trend (for Supplementary Material)
u = sort(dvec(dvec(:,2) > 0, 2)); % Non-zero casualties
th = u(1:end-2); % Thresholds
lambda_pot = zeros(size(th));
for i = 1:length(th)
    iy = find(dvec(:,2) > th(i));
    int = diff(dvec(iy,1));
    lambda_pot(i) = 1 / mean(int); % Rate parameter λ = 1 / mean inter-arrival
end

figure;
plot(th, lambda_pot, 'b-o', 'LineWidth', 1.5);
xlabel('Threshold Casualties (u)');
ylabel('Rate Parameter λ (per year)');
title('Thinning of Rate Parameter λ with Increasing Casualty Threshold (POT)');
grid on;
legend('λ vs. Threshold', 'Location', 'best');
saveas(gcf, 'POTRateParameterTrend.png', 'png');
% Note: Decreasing λ reflects thinned Poisson rate with higher thresholds (Embrechts et al., 1997).
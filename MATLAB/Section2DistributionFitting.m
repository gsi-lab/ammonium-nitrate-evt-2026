%% Severity Analysis: Distribution Fitting
% Streamlined script: Fits GPD, Lognormal, Normal; Monte Carlo uncertainty; visual comparison; hypothesis tests
% Data: 25 non-zero casualties (1916-2022)
% Assumptions: Statistics Toolbox available; no custom MLE yet

close all
clear

% Load data (replace with your load ANdata.mat if needed)
load ANdata.mat
nonzero_casualties = dvec(dvec(:,2) > 0, 2); % 25 non-zero casualties
xraw = nonzero_casualties; % Use raw data since theta=0 is fitted
n = length(xraw);
L = 0; % Fixed location/lower bound
H = 10000; % Upper bound for truncation (Bhopal-scale)
H2 = 8.19e9; % World population bound

% Empirical CDF for plotting
[Fi, Xi] = ecdf(xraw);

%% Fit GPD and Monte Carlo Uncertainty
pd_gp = fitdist(xraw, 'GeneralizedPareto', 'Theta', L);
gp_k = pd_gp.k;
gp_sigma = pd_gp.sigma;
gp_theta = pd_gp.theta;

% Log-likelihood for tests
logL_gp = sum(log(pdf(pd_gp, xraw)));

% Truncated GPD for bounded mean
pd_gp_bounded = truncate(pd_gp, L, H);
pd_gp_bounded2 = truncate(pd_gp, L, H2);
meantruncGPD = mean(pd_gp_bounded); % Shadow mean H=10000
meantruncGPD2 = mean(pd_gp_bounded2); % Shadow mean H=8.19e9

% Monte Carlo for 95% CI on CDF
pcov_gp = pd_gp.ParameterCovariance(1:2,1:2); % Covariance for k and sigma
params_gp = [gp_k, gp_sigma];
Sample_gp = mvnrnd(params_gp, pcov_gp, 10000);
Sample_gp(Sample_gp(:,1) < 0 | Sample_gp(:,2) <= 0, :) = []; % Filter invalid

z = linspace(min(xraw), max(xraw) + 10, 100); % x-range for CDF
ymc_gp = zeros(length(z), size(Sample_gp, 1));
for i = 1:size(Sample_gp, 1)
    pd_temp = makedist('GeneralizedPareto', 'k', Sample_gp(i,1), 'sigma', Sample_gp(i,2), 'theta', L);
    ymc_gp(:,i) = cdf(pd_temp, z);
end

% Plot GPD CDF with 95% CI
fig_gp = figure('Visible', 'on');
plot(Xi, Fi, 'ko');
hold on
plot(z, mean(ymc_gp, 2), 'k', z, prctile(ymc_gp, 2.5, 2), 'r', z, prctile(ymc_gp, 97.5, 2), 'r');
hold off
title('Generalized Pareto with 95% Uncertainty on Predictions');
legend('Data', 'Generalized Pareto', '95% CI', 'Location', 'southeast');
xlabel('x (# of Casualties in AN Accidents)');
ylabel('F_X(x) = P(X < x)');
text(300, 0.45, sprintf('with \\xi=%.2f, \\sigma=%.2f, \\theta=%.2f', gp_k, gp_sigma, gp_theta), 'Interpreter', 'tex');
exportgraphics(fig_gp, 'GPD_fit.png', 'Resolution', 300);

%% Fit Lognormal and Monte Carlo Uncertainty
pd_logn = fitdist(xraw, 'Lognormal');
ln_mu = pd_logn.mu;
ln_sigma = pd_logn.sigma;

% Log-likelihood for tests
logL_logn = sum(log(pdf(pd_logn, xraw)));

% Truncated Lognormal for bounded mean
pd_logn_bounded = truncate(pd_logn, L, H);
meantruncLogn = mean(pd_logn_bounded);

% Monte Carlo for 95% CI on CDF (using parameter covariance)
pcov_ln = pd_logn.ParameterCovariance;
params_ln = [ln_mu, ln_sigma];
Sample_ln = mvnrnd(params_ln, pcov_ln, 10000);
Sample_ln(Sample_ln(:,2) <= 0, :) = []; % Filter invalid sigma

ymc_ln = zeros(length(z), size(Sample_ln, 1));
for i = 1:size(Sample_ln, 1)
    pd_temp = makedist('Lognormal', 'mu', Sample_ln(i,1), 'sigma', Sample_ln(i,2));
    ymc_ln(:,i) = cdf(pd_temp, z);
end

% Plot Lognormal CDF with 95% CI
fig_ln = figure('Visible', 'on');
plot(Xi, Fi, 'ko');
hold on
plot(z, mean(ymc_ln, 2), 'k', z, prctile(ymc_ln, 2.5, 2), 'r', z, prctile(ymc_ln, 97.5, 2), 'r');
hold off
title('Lognormal with 95% Uncertainty on Predictions');
legend('Data', 'Lognormal', '95% CI', 'Location', 'southeast');
xlabel('x (# of Casualties in AN Accidents)');
ylabel('F_X(x) = P(X < x)');
text(300, 0.45, sprintf('with \\mu=%.2f, \\sigma=%.2f', ln_mu, ln_sigma), 'Interpreter', 'tex');
exportgraphics(fig_ln, 'Lognormal_fit.png', 'Resolution', 300);

%% Fit Normal and Monte Carlo Uncertainty
pd_norm = fitdist(xraw, 'Normal');
norm_mu = pd_norm.mu;
norm_sigma = pd_norm.sigma;

% Log-likelihood for tests (optional, for completeness)
logL_norm = sum(log(pdf(pd_norm, xraw)));

% Truncated Normal for bounded mean (though not typically heavy-tailed)
pd_norm_bounded = truncate(pd_norm, L, H);
meantruncNorm = mean(pd_norm_bounded);

% Monte Carlo for 95% CI on CDF (using parameter covariance)
pcov_norm = pd_norm.ParameterCovariance;
params_norm = [norm_mu, norm_sigma];
Sample_norm = mvnrnd(params_norm, pcov_norm, 10000);
Sample_norm(Sample_norm(:,2) <= 0, :) = []; % Filter invalid sigma

ymc_norm = zeros(length(z), size(Sample_norm, 1));
for i = 1:size(Sample_norm, 1)
    pd_temp = makedist('Normal', 'mu', Sample_norm(i,1), 'sigma', Sample_norm(i,2));
    ymc_norm(:,i) = cdf(pd_temp, z);
end

% Plot Normal CDF with 95% CI
fig_norm = figure('Visible', 'on');
plot(Xi, Fi, 'ko');
hold on
plot(z, mean(ymc_norm, 2), 'k', z, prctile(ymc_norm, 2.5, 2), 'r', z, prctile(ymc_norm, 97.5, 2), 'r');
hold off
title('Normal with 95% Uncertainty on Predictions');
legend('Data', 'Normal', '95% CI', 'Location', 'southeast');
xlabel('x (# of Casualties in AN Accidents)');
ylabel('F_X(x) = P(X < x)');
text(300, 0.45, sprintf('with \\mu=%.2f, \\sigma=%.2f', norm_mu, norm_sigma), 'Interpreter', 'tex');
exportgraphics(fig_norm, 'Normal_fit.png', 'Resolution', 300);

%% Visual Comparison Figure with Tiled Layout
fig_compare = figure('Visible', 'on');
t = tiledlayout(3,1);

% Subplot 1: GPD
nexttile
plot(Xi, Fi, 'ko');
hold on
plot(z, mean(ymc_gp, 2), 'k', z, prctile(ymc_gp, 2.5, 2), 'r', z, prctile(ymc_gp, 97.5, 2), 'r');
hold off
title('GPD Fit');
legend('Data', 'Fit', '95% CI', 'Location', 'southeast');
xlabel('x (# of Casualties)');
ylabel('CDF F(x)');

% Subplot 2: Lognormal
nexttile
plot(Xi, Fi, 'ko');
hold on
plot(z, mean(ymc_ln, 2), 'k', z, prctile(ymc_ln, 2.5, 2), 'r', z, prctile(ymc_ln, 97.5, 2), 'r');
hold off
title('Lognormal Fit');
legend('Data', 'Fit', '95% CI', 'Location', 'southeast');
xlabel('x (# of Casualties)');

% Subplot 3: Normal
nexttile
plot(Xi, Fi, 'ko');
hold on
plot(z, mean(ymc_norm, 2), 'k', z, prctile(ymc_norm, 2.5, 2), 'r', z, prctile(ymc_norm, 97.5, 2), 'r');
hold off
title('Normal Fit');
legend('Data', 'Fit', '95% CI', 'Location', 'southeast');
xlabel('x (# of Casualties)');

% Export the comparison figure
exportgraphics(fig_compare, 'CDF_comparison.png', 'Resolution', 300);

%% Hypothesis Testing and Statistics
%% Hypothesis Testing and Statistics
% KS Tests
[h_gp, p_gp] = kstest(xraw, 'CDF', pd_gp);
[h_ln, p_ln] = kstest(xraw, 'CDF', pd_logn);
[h_norm, p_norm] = kstest(xraw, 'CDF', pd_norm);
% Anderson-Darling Tests
[h_ad_gp, p_ad_gp] = adtest(xraw, 'Distribution', pd_gp, 'Alpha', 0.05);
[h_ad_ln, p_ad_ln] = adtest(xraw, 'Distribution', pd_logn, 'Alpha', 0.05);
[h_ad_norm, p_ad_norm] = adtest(xraw, 'Distribution', pd_norm, 'Alpha', 0.05);
% Shapiro-Wilk test
%[h_sw, p_sw, W_sw] = swtest(xraw);
% AIC and BIC for all
k_gp = 2; k_ln = 2; k_norm = 2;
AIC_gp = -2 * logL_gp + 2 * k_gp;
BIC_gp = -2 * logL_gp + k_gp * log(n);
AIC_ln = -2 * logL_logn + 2 * k_ln;
BIC_ln = -2 * logL_logn + k_ln * log(n);
AIC_norm = -2 * logL_norm + 2 * k_norm;
BIC_norm = -2 * logL_norm + k_norm * log(n);
% Vuong Test for GPD vs. Lognormal
d = log(pdf(pd_gp, xraw)) - log(pdf(pd_logn, xraw)); % Log-likelihood difference per observation
d_bar = mean(d); % Mean difference
var_d = var(d); % Variance of differences
V = d_bar / sqrt(var_d / n); % Vuong statistic
p_value_vuong = 2 * (1 - normcdf(abs(V))); % Two-tailed p-value
disp('KS p-values: GPD = ' + string(p_gp) + ', Lognormal = ' + string(p_ln) + ', Normal = ' + string(p_norm));
disp('AD p-values: GPD = ' + string(p_ad_gp) + ', Lognormal = ' + string(p_ad_ln) + ', Normal = ' + string(p_ad_norm));
%disp('Shapiro-Wilk p-value for Normal = ' + string(p_sw));
disp('Vuong Test (GPD vs Lognormal): V = ' + string(V) + ', p-value = ' + string(p_value_vuong));
if abs(V) < 1.96
    disp('GPD and Lognormal are not statistically distinguishable at 5% level.');
else
    if V > 0
        disp('GPD is preferred over Lognormal.');
    else
        disp('Lognormal is preferred over GPD.');
    end
end
disp('AIC: GPD = ' + string(AIC_gp) + ', Lognormal = ' + string(AIC_ln) + ', Normal = ' + string(AIC_norm));
disp('BIC: GPD = ' + string(BIC_gp) + ', Lognormal = ' + string(BIC_ln) + ', Normal = ' + string(BIC_norm));
disp('Truncated Means: GPD (H=10000) = ' + string(meantruncGPD) + ', GPD (H=8.19e9) = ' + string(meantruncGPD2) + ', Lognormal (H=10000) = ' + string(meantruncLogn));

%% Publication-Ready Comparison Figure
% This figure overlays the two best fits (GPD and Lognormal) for direct comparison
% and includes an inset plot to highlight the behavior in the extreme upper tail.

fig_pub = figure('Position', [100, 100, 800, 600]); % Create a nicely sized figure

% --- Main Plot (CDF Comparison) ---
h_main = axes('Position', [0.1 0.1 0.85 0.85]); % Main axes
hold(h_main, 'on');

% Plot Empirical Data
plot(h_main, Xi, Fi, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Empirical Data');

% Plot Fitted Models
plot(h_main, z, mean(ymc_gp, 2), 'b-', 'LineWidth', 2, 'DisplayName', 'GPD Fit');
plot(h_main, z, mean(ymc_ln, 2), 'r--', 'LineWidth', 2, 'DisplayName', 'Lognormal Fit');
plot(h_main, z, mean(ymc_norm, 2), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Normal Fit');

% Formatting
title(h_main, 'Comparison of Distribution Fits for AN Accident Casualties');
xlabel(h_main, 'x (# of Casualties)');
ylabel(h_main, 'Cumulative Probability F(x)');
legend(h_main, 'show', 'Location', 'southeast');
grid(h_main, 'on');
box(h_main, 'on');
set(h_main, 'FontSize', 12);

% --- Insert Plot (Tail Behavior) ---
% This shows the survival function (1-CDF) on a log scale to emphasize the tail

fig_pub = figure('Position', [100, 100, 800, 600]); % Create a nicely sized figure

% --- Main Plot (CDF Comparison) ---
h_main = axes('Position', [0.1 0.1 0.85 0.85]);
hold(h_main, 'on');

% Plot Empirical Data and Fitted Models (same as before)
plot(h_main, Xi, Fi, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Empirical Data');
plot(h_main, z, mean(ymc_gp, 2), 'b-', 'LineWidth', 2, 'DisplayName', 'GPD Fit');
plot(h_main, z, mean(ymc_ln, 2), 'r--', 'LineWidth', 2, 'DisplayName', 'Lognormal Fit');
plot(h_main, z, mean(ymc_norm, 2), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', 'Normal Fit');

% Formatting (same as before)
title(h_main, 'Comparison of Distribution Fits for AN Accident Casualties');
xlabel(h_main, 'x (# of Casualties)');
ylabel(h_main, 'Cumulative Probability F(x)');
legend(h_main, 'show', 'Location', 'Southwest');
legend('boxoff')
grid(h_main, 'on');
box(h_main, 'on');
set(h_main, 'FontSize', 12);

% --- Inset Plot (Tail Behavior) ---
% Positioned in the lower-right to be non-obstructive
h_inset = axes('Position', [0.55, 0.2, 0.35, 0.35]); % <-- CHANGED LINE: Adjusted position
hold(h_inset, 'on');

% Plotting for inset (same as before)
semilogy(h_inset, Xi, 1-Fi, 'ko', 'MarkerFaceColor', 'k');
semilogy(h_inset, z, 1-mean(ymc_gp, 2), 'b-', 'LineWidth', 2);
semilogy(h_inset, z, 1-mean(ymc_ln, 2), 'r--', 'LineWidth', 2);
semilogy(h_inset, z, 1-mean(ymc_norm, 2), 'k:', 'LineWidth', 2);

% Formatting Inset (same as before)
title(h_inset, 'Tail Behavior (Log Scale)');
xlabel(h_inset, 'x');
ylabel(h_inset, 'Survival F. 1-F(x)');
grid(h_inset, 'on');
box(h_inset, 'on');
xlim(h_inset, [quantile(xraw, 0.5), max(xraw) + 50]);
set(h_inset, 'FontSize', 10);

hold(h_main, 'off');
hold(h_inset, 'off');

% --- Export the Final Figure ---
exportgraphics(fig_pub, 'CDF_Comparison_Publication.png', 'Resolution', 300);
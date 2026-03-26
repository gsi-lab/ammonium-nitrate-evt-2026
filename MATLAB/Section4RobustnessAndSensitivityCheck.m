% Robustness and Sensitivity Analysis for AN Accident Casualties
% Implements bootstrap resampling and Monte Carlo sensitivity for GPD shape parameter xi
% Data: 25 non-zero casualty events from AN accidents (1916-2022)
% Outputs: Bootstrap xi distribution (mean, mode, 95% CI), Monte Carlo xi stability

clear all; close all; clc;
load('ANdata.mat'); % Assumes ANdata.mat is in the path
all_events = dvec; % Contains [DecimalYear, Casualties] for all 35 events
casualties_nonzero = all_events(all_events(:,2) > 0, 2); % n=25 non-zero 
% Input data: 25 non-zero casualty events from Supplementary Table S1
data = casualties_nonzero;
n = length(data); % Sample size: 25

% Bootstrap settings
B = 10000; % Number of bootstrap resamples
xi_boot = zeros(B, 1); % Store bootstrap xi estimates

% Bootstrap resampling for GPD shape parameter xi
rng(123); % Set seed for reproducibility
for i = 1:B
    % Resample with replacement
    boot_sample = randsample(data, n, true);
    % Fit GPD using MLE
    pd = fitdist(boot_sample, 'GeneralizedPareto', 'Theta', 0);
    xi_boot(i) = pd.k; % Shape parameter xi
end

% Calculate bootstrap statistics
xi_mean = mean(xi_boot);
xi_std = std(xi_boot);
xi_ci = prctile(xi_boot, [2.5, 97.5]); % 95% confidence interval

% Estimate mode using kernel density
[xi_density, xi_grid] = ksdensity(xi_boot, 'Support', [min(xi_boot)-0.5, max(xi_boot)+0.5], 'NumPoints', 1000);
[~, idx] = max(xi_density);
xi_mode = xi_grid(idx);

% Monte Carlo sensitivity analysis (omit up to 20% of data)
max_omit = floor(0.2 * n); % Max 5 observations omitted
n_mc = 10000; % Monte Carlo iterations
xi_mc = zeros(n_mc, max_omit); % Store xi for each omission level
omit_counts = 1:max_omit; % Number of observations to omit (1 to 5)

for j = 1:max_omit
    for i = 1:n_mc
        % Randomly omit j observations
        idx = randperm(n, n-j); % Indices to keep
        mc_sample = data(idx);
        % Fit GPD to reduced sample
        pd_mc = fitdist(mc_sample, 'GeneralizedPareto', 'Theta', 0);
        xi_mc(i, j) = pd_mc.k;
    end
end

% Calculate Monte Carlo statistics
xi_mc_mean = mean(xi_mc);
xi_mc_std = std(xi_mc);
xi_mc_ci = prctile(xi_mc, [2.5, 97.5], 1); % 95% CI for each omission level

% Display results
fprintf('Bootstrap Results for xi (n=%d, B=%d):\n', n, B);
fprintf('Mean xi: %.4f\n', xi_mean);
fprintf('Std xi: %.4f\n', xi_std);
fprintf('95%% CI: [%.4f, %.4f]\n', xi_ci(1), xi_ci(2));
fprintf('Mode xi: %.4f\n', xi_mode);
fprintf('\nMonte Carlo Sensitivity (omit 1 to %d observations):\n', max_omit);
for j = 1:max_omit
    fprintf('Omit %d: Mean xi = %.4f, Std xi = %.4f, 95%% CI = [%.4f, %.4f]\n', ...
        j, xi_mc_mean(j), xi_mc_std(j), xi_mc_ci(1,j), xi_mc_ci(2,j));
end

% Plot bootstrap xi distribution
figure;
histogram(xi_boot, 50, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'k');
hold on;
plot(xi_grid, xi_density, 'r-', 'LineWidth', 2);
xline(xi_mode, 'g--', sprintf('Mode xi = %.2f', xi_mode), 'LineWidth', 2);
xline(xi_mean, 'k--', sprintf('Mean xi = %.2f', xi_mean), 'LineWidth', 2);
title('Bootstrap Distribution of GPD Shape Parameter \xi');
xlabel('\xi'); ylabel('Density');
legend('Histogram', 'Kernel Density', 'Mode', 'Mean');
grid on;
exportgraphics(gcf, 'BootsrapdistributionGPDXi.png', 'Resolution', 300);

% Plot Monte Carlo xi stability
figure;
errorbar(omit_counts, xi_mc_mean, xi_mc_std, 'o-', 'LineWidth', 2);
title('Monte Carlo Sensitivity: GPD \xi vs. Number of Omitted Observations');
xlabel('Number of Omitted Observations'); ylabel('Mean \xi ± Std');
grid on;
exportgraphics(gcf, 'MonteCarloSensitivity.png', 'Resolution', 300);
% Save results

%% Improved Visualization with Tiledlayout for AN Robustness Results
% Requires xi_boot, xi_mc, xi_mean, xi_mode, xi_ci, xi_mc_mean, xi_mc_std from previous run
% Uses exportgraphics for high-quality output
figure
% Create tiled layout with increased spacing
t = tiledlayout(1, 2, 'TileSpacing', 'normal', 'Padding', 'normal');

% Tile 1: Bootstrap Distribution of xi (Enhanced)
%nexttile([1 1.5]); % Allocate more width to bootstrap plot
nexttile
histogram(xi_boot, 75, 'Normalization', 'pdf', 'FaceColor', [0.5 0.5 1], 'EdgeColor', 'k', 'FaceAlpha', 0.6);
hold on;
%[xi_density, xi_grid] = ksdensity(xi_boot, 'Support', [-0.5, 2.5], 'NumPoints', 1000);
plot(xi_grid, xi_density, 'r-', 'LineWidth', 2);
xline(xi_mode, 'g--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
text(xi_mode + 0.1, 0.8, sprintf('Mode = %.2f', xi_mode), 'Color', 'g', 'FontSize', 10);
xline(xi_mean, 'k--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
text(xi_mean + 0.1, 0.7, sprintf('Mean = %.2f', xi_mean), 'Color', 'k', 'FontSize', 10);
xline(xi_ci(1), 'b--', 'LineWidth', 1);
xline(xi_ci(2), 'b--', 'LineWidth', 1);
text(xi_ci(1) - 0.2, 0.6, sprintf('CI: [%.2f', xi_ci(1)), 'Color', 'b', 'FontSize', 10);
text(xi_ci(2) + 0.1, 0.6, sprintf('%.2f]', xi_ci(2)), 'Color', 'b', 'FontSize', 10);
title('Bootstrap Distribution of GPD Shape Parameter \xi', 'FontSize', 11);
xlabel('\xi', 'FontSize', 10); ylabel('Density', 'FontSize', 10);
legend('Histogram', 'Kernel Density', 'Mode', 'Mean', '95% CI', 'Location', 'northeastoutside', 'FontSize', 8);
grid on;
xlim([-0.5 2.5]); ylim([0 1.5]);
legend('boxoff')
% Tile 2: Monte Carlo Sensitivity (Enhanced)
nexttile;
errorbar(1:5, xi_mc_mean, xi_mc_std, 'o-', 'LineWidth', 2, 'Color', [0 0.5 0], 'MarkerSize', 8, 'MarkerFaceColor', [0 0.5 0]);
hold on;
plot(1:5, xi_mc_mean, 'b-', 'LineWidth', 1.5);
title('Monte Carlo Sensitivity: GPD \xi vs. Omitted Observations', 'FontSize', 11);
xlabel('Number of Omitted Observations', 'FontSize', 10); ylabel('Mean \xi ± Std', 'FontSize', 10);
xlim([0.5 5.5]); ylim([0.9 1.4]);
grid on;
for i = 1:5
    text(i, xi_mc_mean(i) + 0.02, sprintf('%.2f', xi_mc_mean(i)), 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 10);
end

% Adjust layout
%title(t, 'Robustness and Sensitivity Analysis of AN Accident Data', 'FontSize', 14);
set(gcf, 'Position', [100 100 1000 400]); % Increased figure size

% Export figure with exportgraphics
exportgraphics(gcf, 'AN_robustness_sensitivity_tiled.png', 'Resolution', 300);


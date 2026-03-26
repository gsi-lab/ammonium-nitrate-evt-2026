% Script: AN Accident Risk Analysis - Final T/N Monitoring Plot (50-yr Horizon)
clear all; close all;

% 1. Parameters from your EVT analysis 
lambda = 0.229;     % Overall hazard rate (events per year)
xi = 1.13;          % GPD Shape Parameter
sigma = 24.09;      % GPD Scale Parameter
theta = 0;          % Location threshold

% 2. Monitoring Range
N_ext = logspace(0, 4, 200); 

% 3. Calculating the Systemic Health Monitor 
% S(N) = (1 + xi*N/sigma)^(-1/xi)
S_N_ext = (1 + xi * (N_ext - theta) / sigma).^(-1/xi); 
% T(N) = 1 / (lambda * S(N))
T_N_ext = 1 ./ (lambda * S_N_ext);

% 4. Visualization
figure('Color', 'w', 'Position', [100, 100, 800, 550]);
loglog(N_ext, T_N_ext, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Systemic Health Monitor');
hold on;

% Add Monitoring Horizons
yline(50, 'g--', 'LineWidth', 1.5, 'DisplayName', '50-yr Horizon'); % Updated to 50-yr
yline(100, 'r--', 'LineWidth', 1.5, 'DisplayName', '100-yr Horizon');

% Find 50-year Return Level intersection
[~, idx50] = min(abs(T_N_ext - 50));
plot(N_ext(idx50), 50, 'go', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(N_ext(idx50), 65, sprintf(' N approx %.0f', N_ext(idx50)), ...
    'Color', 'g', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Find 100-year Return Level intersection
[~, idx100] = min(abs(T_N_ext - 100));
plot(N_ext(idx100), 100, 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(N_ext(idx100), 140, sprintf(' N approx %.0f', N_ext(idx100)), ...
    'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Labels and Grid
xlabel('Fatalities (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Return Period T(N) (years)', 'FontSize', 12, 'FontWeight', 'bold');
title('Monitoring Systemic Safety: T/N Return Levels', 'FontSize', 14);
grid on; set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');

% Legend
legend('Location', 'northwest', 'Box', 'on', 'EdgeColor', 'k', 'FontSize', 10);

% 5. Equation Text Box (Lower-Right Position)
str = {'Survival Function: S(N) = (1 + \xi N / \sigma)^{-1/\xi}', ...
       'Return Period: T(N) = 1 / (\lambda S(N))', ...
       'Baseline: \xi = 1.13,\sigma = 24.09, \lambda = 0.23'};
annotation('textbox', [0.55, 0.15, 0.33, 0.15], 'String', str, ...
           'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FaceAlpha', 0.9);

saveas(gcf, 'Systemic_Monitoring_TN_50yr.png');
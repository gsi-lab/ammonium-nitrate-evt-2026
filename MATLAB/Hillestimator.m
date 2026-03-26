%% running the script
% clear 
% close all

load ANdata.mat
nonzero_casualties = dvec(dvec(:,2) > 0, 2); % 25 non-zero casualties
xraw = nonzero_casualties; % Use raw data since
Xo=xraw;
Xj=sort(Xo,'descend'); % sorted data
rj=1:length(Xj); %rank of the data
n=length(Xj);

%% Hill estimator for different u
k=3:n-1; % ignore the largest observation
Hill=[];Hillsd=[];
k=round(k);
for j=1:length(k) % select different kmax that is a certain sample of the tail events.
    u(j)=Xj(k(j)+1);
    ik=k(j);
    [Hill(j) Hsd(j)]=TailHill(Xj,ik); % this gives alpha, tail index
end
mdl = fitlm(k,Hill,'Weights',k,'Intercept',true)%
alfa= mdl.Coefficients.Estimate(1) ; % use weighted approach to consider the small sample size


% Target the third subplot (tile 3)
ax1=nexttile(1);

% Primary axis (ax1) for Hill estimator vs. k in the third subplot
plot(ax1, k, Hill, 'b', k, Hill + 1.96*Hsd, 'r-', k, Hill - 1.96*Hsd, 'r-');
hold(ax1, 'on');
% yline(ax1, alfa, 'k--', 'LineWidth', 1, 'Label', ['\alpha_{weighted}=', num2str(round(alfa,2))]);
xlim(ax1, [0 k(end)]);
ylim(ax1, [0 4]);
xt = k; % k values as x-ticks
xticks(ax1, xt);
ylabel(ax1, 'Hill estimator \alpha');
xlabel(ax1, 'k order statistics');
grid(ax1, 'off');
axis(ax1, 'tight');
legend(ax1, {'\alpha Hill estimator', '\pm 95 CI of Hill estimator'}, 'box', 'off');
hold on
% Secondary axis (ax2) for threshold u vs. k in the same tile
ax2 = axes(t); % Create secondary axes in the same layout
h = plot(ax2, k, Hill, '-'); % Plot same data, but will be hidden
h.Color = 'none'; % Hide the line, use ax2 for labels only
ax2.XAxisLocation = 'top'; % Move x-axis to top
ax2.YAxisLocation = 'right'; % Move y-axis to right (optional, can hide)
ax2.Color = 'none'; % Make background transparent
ax1.Box = 'off'; % Turn off box on primary to avoid overlap
ax2.Box = 'off'; % Turn off box on secondary
ax2.XLim = ax1.XLim; % Sync x-limits with ax1
ax2.YLim = [0 4]; % Match y-limits for alignment
xticks(ax2, xt); % Use same k ticks
% Map k to thresholds (assuming X is sorted data, e.g., casualties)
X = sort(nonzero_casualties, 'descend'); % Replace with your data vector
xticklabels(ax2, compose('%0.1f', round(X(xt), 1))); % Label with u values
xlabel(ax2, 'threshold u');
axis(ax2, 'tight');

% Save the figure (only the third subplot is affected)
saveas(1, 'GraphicalDiagnosticsSeverity', 'png');
%% data uncertainty: bootstrap
% tail index estimation using Hill estimator
dec=1;
if dec==1

nb=1e4;
Xj=[];Hillb=[];
u=10; % threshold 

for i=1:nb

    Xj =datasample(Xo,n); % random sample with replacement
    Xsorted=sort(Xj,'descend');
    kth=max(find(Xsorted>u)); % find the k upper order for this threshold
    Hillb(i,1)=TailHill(Xj,kth-1); % this gives alpha, tail index

    pd_gp = fitdist(Xsorted, 'GeneralizedPareto');

end
pd_gp = fitdist(Xo, 'GeneralizedPareto');
xi_hat=pd_gp.k;

figure
tiledlayout(1,2)
nexttile
histogram(Hillb)
xlabel('Hill estimator for \alpha')
ylabel('Frequency');
xline(1/xi_hat, 'r-', 'LineWidth', 2,'label',['\alpha=',num2str(1/xi_hat)]);
legend('Tail index \alpha','Reference estimate')
legend('Location','best')
legend('boxoff')
nexttile
histogram(1./Hillb)
xlabel('Shape parameter, \xi = \alpha^{-1}')
ylabel('Frequency');
xline(xi_hat, 'r-', 'LineWidth', 2,'label',['\xi=',num2str(xi_hat)]);
legend('shape parameter \xi','Reference estimate')
legend('Location','best')
legend('boxoff')
saveas(gcf,'BootstrapEstimateAlfa','png')
end




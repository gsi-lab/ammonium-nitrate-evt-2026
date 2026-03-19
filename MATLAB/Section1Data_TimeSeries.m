%% Ammonium Nitrate Accident Analysis - Temporal Section
clear
close all
rng(1976); % Set random seed for reproducibility

%% Load and Prepare Data
T = readtable('MajorAccidents.xlsx', 'Sheet', 'Nitrate_2', 'Range', 'D1:H40'); % 39 rows, excluding header
T.Date = datetime(T.Date, 'InputFormat', 'MMMM d, uuuu'); % Parse full dates
Date=datetime(T.Date,InputFormat="MMMM d, uuuu");
ix=find(T.Deaths>30);

figure(1)
t=tiledlayout(1,1);
ax1 = axes(t);
bar(Date,T.Deaths)
xlim([Date(1) '01-Jan-2030'])
xlabel('Date')
axis normal
ylabel('Casualties')
%xt=linspace(Date(1),Date(end),30);xticks(ax1,xt)
title('Historical records of major accidents involving Ammonium nitrate')
text(Date(ix),T.Deaths(ix),T.Location(ix))
%datetick('x','yyyy')

saveas(1,'TimeSeries','png')

figure(2)
t=tiledlayout(1,1);
ax1 = axes(t);
semilogy(Date,log(T.Deaths),'ko')
xlim([Date(1) '01-Jan-2030'])
xlabel('Date')
axis normal
ylabel('Casualties (log scale)')
%xt=linspace(Date(1),Date(end),30);xticks(ax1,xt)
title('Timeline of Major Ammnonium Nitrate (AN) accidents: 1916-2022')
text(Date(ix),T.Deaths(ix),T.Location(ix))
%datetick('x','yyyy')

saveas(2,'TimeSeriesLog','png')
%% Setup
% Let's start clean :)
clear, clc, close all

% Add dependencies to the path
addpath(genpath('deps'))

% Load data
load('data/2017-09-04-master_dataset.mat')

figPath = 'figs/corr_graphs';
mkdirto(figPath);

%% Setup data
mice = behavior.MouseID(behavior.Region == "WT Control");
raw_metrics = metrics(~(contains(metrics,'_PC') | contains(metrics,'GR_')));
raw = behavior(mice, raw_metrics);

%% WT
% Z-score
mu = nanmean(raw.Metrics);
sigma = nanstd(raw.Metrics);
Z = raw; Z.Metrics = (Z.Metrics - mu) ./ sigma;

% Plot
alpha = 0.05;
plotcorrgraph(Z, alpha)
titlef('WT (n = %d, \\alpha = %.2f)', height(Z), alpha)
export_fig(ff(figPath,'_WT'), '-png', '-eps')

%% Experimentals
alpha = 0.05;
for r = 1:numel(regions)
    region = regions{r};
    for a = 1:numel(ages)
        age = ages{a};
        
        % Pull out mice
        idx = behavior.Region(X.MouseID) == string(region) & behavior.Age(X.MouseID) == string(age);
        Z = X(idx,:);
        
        % Exclude PCs
        Z = Z(:,~contains(Z.Properties.VariableNames,'_PC'));
        
        % Plot
        plotcorrgraph(Z, alpha)
        titlef('%s %s (n = %d, \\alpha = %.2f)', region, age, height(Z), alpha)
        export_fig(ff(figPath,sprintf('%s_%s',region,age)), '-png', '-eps')
        close(gcf)
    end
end

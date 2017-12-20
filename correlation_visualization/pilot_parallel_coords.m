cc
%% Setup
% Load data
load('data/2017-09-04-master_dataset.mat')

%% Setup data
mice = behavior.MouseID(behavior.Region == "WT Control");
raw_metrics = metrics(~(contains(metrics,'_PC') | contains(metrics,'GR_')));
raw = behavior(mice, raw_metrics);

%% WT
% Z-score
mu = nanmean(raw.Metrics);
sigma = nanstd(raw.Metrics);
Z_wt = raw; Z_wt.Metrics = (Z_wt.Metrics - mu) ./ sigma;
[rho_wt,p_wt] = corr(Z_wt.Metrics,'type','Spearman','rows','pairwise');

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
        Z = Z(:,~(contains(Z.Properties.VariableNames,'_PC') | contains(Z.Properties.VariableNames,'GR_')));
        
        [rho,p] = corr(Z.Metrics,'type','Spearman','rows','pairwise');
        
        % Plot
%         plotcorrgraph(Z, alpha, [], [],[], rho_wt)
%         titlef('%s %s (n = %d, \\alpha = %.2f)', region, age, height(Z), alpha)
        return
%         export_fig(ff(figPath,sprintf('%s_%s',region,age)), '-png', '-eps')
%         close(gcf)
    end
end

%%
figure,figclosekey
hold on
% plot(Z_wt{:,1:4}','b-')
% plot(Z{:,1:4}','r-')

plot(Z_wt.Metrics','b-')
plot(Z.Metrics','r-')


% plot(nanmean(Z.Metrics) - nanmean(Z_wt.Metrics))

%%
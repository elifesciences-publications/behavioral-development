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
[Zs,rhos,ps] = replicate(cell(numel(regions),numel(ages)));
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
        
        Zs{r,a} = Z.Metrics;
        rhos{r,a} = rho;
        ps{r,a} = p;
    end
end

%%
labels = Z_wt.Properties.VariableNames;
Z_wt = Z_wt.Metrics;
save('test_corr_data.mat','Z_wt','rho_wt','p_wt','Zs','rhos','ps','regions','ages','labels')

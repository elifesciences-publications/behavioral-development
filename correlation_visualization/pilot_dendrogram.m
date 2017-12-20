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
delta = rho - rho_wt;

Z = linkage(delta);
% T = cluster(Z,'maxclust',5);

figure,figclosekey
subplot(2,1,1)
[h,T,outperm] = dendrogram(Z,0);
labels = Z_wt.Properties.VariableNames(outperm);

subplot(2,1,2)
imagesc(delta(outperm,outperm))
xticks(1:numel(labels)),yticks(1:numel(labels))
xticklabels(labels),xtickangle(90)
yticklabels(labels)
% xticklabels(string(outperm(xticks())))
% yticklabels(string(outperm(yticks())))


%%
% Create graph
G = graph(rho,'upper','OmitSelfLoops');

% Remove missing
G = G.rmedge(find(isnan(G.Edges.Weight)));

% Pull out edge data
w0 = G.Edges.Weight;
pG0 = p(sub2ind(size(p),G.Edges.EndNodes(:,1),G.Edges.EndNodes(:,2)));

% Threshold
% G = G.rmedge(find(pG0 > alpha));
G = G.rmedge(find(pG0 > 0.05));

% Update edge data
w = G.Edges.Weight;
pG = p(sub2ind(size(p),G.Edges.EndNodes(:,1),G.Edges.EndNodes(:,2)));

positive_color = [1 0 0];
negative_color = [0 0 1];

%%
labels = Z.Properties.VariableNames;
[grps,names] = findgroups(extractBefore(labels,'_'));

yd = cellcat(af(@(i)(1:sum(grps==i)) - sum(grps==i)/2, unique(grps)),2);
xd = cellcat(af(@(i)ones(1,sum(grps==i))*i, unique(grps)),2);

figure,figclosekey
plot(G,'XData',xd,'YData',yd,...
    'linewidth',norm2unit(abs(w)) .* 5 + 1,... % width proportional to rho
    'edgecolor',positive_color .* (w < 0) + negative_color .* (w > 0),... % color = sign
    'nodelabel',Z.Properties.VariableNames,... % node names = metrics
    'edgelabel', af(@(x)num2str(x,'%.3f'), w))


%%
labels = Z.Properties.VariableNames;
[grps,names] = findgroups(extractBefore(labels,'_'));
N = arrayfun(@(i)sum(grps==i),unique(grps));

%%
theta = af(@(i)linspace(0,2*pi,N(i)+1),unique(grps));
xd = af(@(i)cos(theta{i}(1:end-1)),unique(grps));
yd = af(@(i)sin(theta{i}(1:end-1)),unique(grps));

xd{1} = xd{1} - 3; yd{1} = yd{1} + 0;
xd{2} = xd{2} + 3; yd{2} = yd{2} + 0;
xd{3} = xd{3} + 0; yd{3} = yd{3} + 3;

xd = cellcat(xd,2); yd = cellcat(yd,2);


figure,figclosekey
plot(G,'XData',xd,'YData',yd,...
    'linewidth',norm2unit(abs(w)) .* 5 + 1,... % width proportional to rho
    'edgecolor',positive_color .* (w < 0) + negative_color .* (w > 0),... % color = sign
    'nodelabel',Z.Properties.VariableNames,... % node names = metrics
    'edgelabel', af(@(x)num2str(x,'%.3f'), w))
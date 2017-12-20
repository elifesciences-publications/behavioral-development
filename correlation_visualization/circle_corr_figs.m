% Goal is to make circle correlation graphs

cc
%% Setup
% Load data
load('data/2017-09-04-master_dataset.mat')

%%
isMotor = {
    'EPM_AntiHesitation'      false
    'EPM_Distance'            true
    'EPM_ExplorationEntrances' false
    'EPM_ExplorationTime'     false
    'EPM_OpenArmPref'         false
    'EPM_Velocity'            true
    'SC_BSDistance'           true
    'SC_NoveltySeeking'       false
    'SC_SocialPreference'     false
    'SC_TestDistance'         true
    'YM_AcqAbility'           false
    'YM_AcqInitialLR'         false
    'YM_AcqSecondaryLR'       false
    'YM_Distance'             true
    'YM_EarlyRevAbility'      false
    'YM_EarlyRevInitialLR'    false
    'YM_EarlyRevSecondaryLR'  false
    'YM_LateRevAbility'       false
    'YM_LateRevInitialLR'     false
    'YM_Velocity'             true
    };

motor_labels = isMotor([isMotor{:,2}],1);
nonmotor_labels = isMotor(~[isMotor{:,2}],1);

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


name = 'WT';

% Save correlations
labels = Z_wt.Properties.VariableNames;
tablePath = ff('data',[name '.xlsx']);
if exists(tablePath); delete(tablePath); end
writetable(Z_wt,tablePath,'WriteRowNames',true,'Sheet','Z')
writetable(array2table(rho_wt,'RowNames',labels,'VariableNames',labels),...
    tablePath,'WriteRowNames',true,'Sheet','rho')
writetable(array2table(p_wt,'RowNames',labels,'VariableNames',labels),...
    tablePath,'WriteRowNames',true,'Sheet','p')
RemoveSheet123(tablePath)

% Plot
% plot_graph(name,Z_wt,rho_wt,p_wt,motor_labels)
% title(replace(name,'_',' / '))
% 
% export_fig(ff('figs',name), '-png', '-eps')
% close(gcf)


%% Experimentals
alpha = 0.05;
for r = 1:numel(regions)
    region = regions{r};
    for a = 1:numel(ages)
        age = ages{a};
        name = [age '_' region];
        
        % Pull out mice
        idx = behavior.Region(X.MouseID) == string(region) & behavior.Age(X.MouseID) == string(age);
        Z = X(idx,:);
        
        % Exclude PCs
        Z = Z(:,~(contains(Z.Properties.VariableNames,'_PC') | contains(Z.Properties.VariableNames,'GR_')));
        
        % Compute correlations
        [rho,p] = corr(Z.Metrics,'type','Spearman','rows','pairwise');
        
        % Save correlations
        labels = Z.Properties.VariableNames;
        tablePath = ff('data',[name '.xlsx']);
        if exists(tablePath); delete(tablePath); end
        writetable(Z,tablePath,'WriteRowNames',true,'Sheet','Z')
        writetable(array2table(rho,'RowNames',labels,'VariableNames',labels),...
            tablePath,'WriteRowNames',true,'Sheet','rho')
        writetable(array2table(p,'RowNames',labels,'VariableNames',labels),...
            tablePath,'WriteRowNames',true,'Sheet','p')
        RemoveSheet123(tablePath)
        
        % Plot
%         plot_graph(name,Z,rho,p,motor_labels)
%         title(replace(name,'_',' / '))
%         
%         export_fig(ff('figs',sprintf('%s_%s',region,age)), '-png', '-eps')
%         close(gcf)
    end
end


%% 
function plot_graph(name,Z,rho,p,motor_labels)
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


labels = Z.Properties.VariableNames;
[grps,names] = findgroups(extractBefore(labels,'_'));
N = arrayfun(@(i)sum(grps==i),unique(grps));

%
theta = af(@(i)linspace(0,2*pi,N(i)+1),unique(grps));
xd = af(@(i)cos(theta{i}(1:end-1)),unique(grps));
yd = af(@(i)sin(theta{i}(1:end-1)),unique(grps));

xd{1} = xd{1} - 3; yd{1} = yd{1} + 0;
xd{2} = xd{2} + 3; yd{2} = yd{2} + 0;
xd{3} = xd{3} + 0; yd{3} = yd{3} + 3;

xd = cellcat(xd,2); yd = cellcat(yd,2);


figure,figclosekey
figsize(1000,750)
h = plot(G,'XData',xd,'YData',yd,...
    'linewidth',norm2unit(abs(w)) .* 5 + 1,... % width proportional to rho
    'edgecolor',positive_color .* (w < 0) + negative_color .* (w > 0),... % color = sign
    'nodelabel',1:size(rho,1),... % node names = numbers
    'edgelabel', af(@(x)num2str(x,'%.2f'), w));
%     'nodelabel',Z.Properties.VariableNames,... % node names = metrics
%     'NodeCData',contains(labels,motor_labels)+1,...
% colormap(viridis())

cm = colormap(viridis(2));
% highlight(h,find(contains(labels,motor_labels)),'Marker','o','MarkerSize',15,'NodeColor',cm(1,:))
% highlight(h,find(~contains(labels,motor_labels)),'Marker','o','MarkerSize',15,'NodeColor',cm(2,:))
% highlight(h,find(~contains(labels,motor_labels)),'Marker','o','MarkerSize',15,'NodeColor',[1 1 1])

hold on
hnodes = gobjects(2,1);
hnodes(1) = plot(xd(contains(labels,motor_labels)),yd(contains(labels,motor_labels)),'o','markersize',12,'markerfacecolor',[1 1 1],'markeredgecolor',cm(1,:));
hnodes(2) = plot(xd(~contains(labels,motor_labels)),yd(~contains(labels,motor_labels)),'o','markersize',12,'markerfacecolor',cm(2,:),'markeredgecolor',cm(1,:));

legend(hnodes,{'Motor','Non-motor'})
fontsize(16)

text(0,1,af(@(i)sprintf('%d. %s',i,labels{i}),1:numel(labels)),'units','normalized','VerticalAlignment','top','FontSize',8,'interpreter','none')

painters
noax

end
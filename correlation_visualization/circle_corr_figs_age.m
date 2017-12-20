% Goal here is to make circle corr figs with sig diffs
cc
%% Setup
% Load data
load('data/2017-09-04-master_dataset.mat')

% Metrics
raw_metrics = metrics(~(contains(metrics,'_PC') | contains(metrics,'GR_')));

% Figures
figDir = 'figs/v4_age'; mkdirto(figDir)

%% WT
mice = behavior.MouseID(behavior.Region == "WT Control");
raw = behavior(mice, raw_metrics);

% Z-score
mu = nanmean(raw.Metrics);
sigma = nanstd(raw.Metrics);
Z = raw; Z.Metrics = (Z.Metrics - mu) ./ sigma;
[wt_rho,wt_p] = corr(Z.Metrics,'type','Spearman','rows','pairwise');

plot_corrs(wt_rho,wt_p,raw_metrics)
title('WT','FontSize',18)

export_fig(ff(figDir,'WT'),'-eps','-png'); close(gcf)

%% Experimentals
for a = 1:numel(ages)
    age = ages{a};
    name = age;

    % Pull out mice
    idx = behavior.Age(X.MouseID) == string(age);
    Z = X(idx,:);

    % Exclude PCs
    Z = Z(:,~(contains(Z.Properties.VariableNames,'_PC') | contains(Z.Properties.VariableNames,'GR_')));

    % Compute correlations
    [rho,p] = corr(Z.Metrics,'type','Spearman','rows','pairwise');

    plot_corrs(rho,p,raw_metrics)
    title(name,'FontSize',18)

    export_fig(ff(figDir,name),'-eps','-png'); close(gcf)
end

%% Experimentals (diff)
alpha = 0.05;
p_mask = wt_p < alpha;
figDir_diff = ff(figDir,'diff'); mkdirto(figDir_diff);

for a = 1:numel(ages)
    age = ages{a};
    name = [age '_diff'];

    % Pull out mice
    idx = behavior.Age(X.MouseID) == string(age);
    Z = X(idx,:);

    % Exclude PCs
    Z = Z(:,~(contains(Z.Properties.VariableNames,'_PC') | contains(Z.Properties.VariableNames,'GR_')));

    % Compute correlations
    [rho,p] = corr(Z.Metrics,'type','Spearman','rows','pairwise');

    % Mask
    p(p_mask) = 1;

    % Plot
    plot_corrs(rho,p,raw_metrics)
    title(sprintf('%s (n.s. in WT)', age),'FontSize',18)

    % Save
    export_fig(ff(figDir_diff,name),'-eps','-png'); close(gcf)
end


%% Helpers
function plot_corrs(rho,p,labels)
%%% Node locations
node_locs = NaN(numel(labels),2);
node_cols = zeros(numel(labels),3);
cm = viridis(3) .* 0.7;

assay = 'EPM';
ctr = [-3 0];
idx = contains(labels,assay);
dtheta = 2*pi / sum(idx);
theta = 0:dtheta:(dtheta*(sum(idx)-1));
node_locs(idx,:) = [cos(theta(:)), sin(theta(:))] + ctr;
node_cols(idx,:) = node_cols(idx,:) + cm(1,:);

assay = 'SC';
ctr = [3 0];
idx = contains(labels,assay);
dtheta = 2*pi / sum(idx);
theta = 0:dtheta:(dtheta*(sum(idx)-1)); theta = theta + deg2rad(20);
node_locs(idx,:) = [cos(theta(:)), sin(theta(:))] + ctr;
node_cols(idx,:) = node_cols(idx,:) + cm(2,:);

assay = 'YM';
ctr = [0 3];
idx = contains(labels,assay);
dtheta = 2*pi / sum(idx);
theta = 0:dtheta:(dtheta*(sum(idx)-1));
node_locs(idx,:) = [cos(theta(:)), sin(theta(:))] + ctr;
node_cols(idx,:) = node_cols(idx,:) + cm(3,:);


% Motor/nonmotor metric categories
motor_metrics = { 
    'EPM_Distance'   
    'EPM_Velocity'   
    'SC_BSDistance'  
    'SC_TestDistance'
    'YM_Distance'    
    'YM_Velocity'    };
% nonmotor_metrics = {
%     'EPM_AntiHesitation'      
%     'EPM_ExplorationEntrances'
%     'EPM_ExplorationTime'     
%     'EPM_OpenArmPref'         
%     'SC_NoveltySeeking'       
%     'SC_SocialPreference'     
%     'YM_AcqAbility'           
%     'YM_AcqInitialLR'         
%     'YM_AcqSecondaryLR'       
%     'YM_EarlyRevAbility'      
%     'YM_EarlyRevInitialLR'    
%     'YM_EarlyRevSecondaryLR'  
%     'YM_LateRevAbility'       
%     'YM_LateRevInitialLR'     };

isMotor = contains(labels,motor_metrics);
% node_cols = zeros(numel(labels),3);
% node_cols(isMotor,:) = node_cols(isMotor,:) + motor_col;
% node_cols(~isMotor,:) = node_cols(~isMotor,:) + nonmotor_col;

%%% Plot
figure,figclosekey
figsize(800,650)
hold on

% Edges
edge_scale = 5;
alpha = 0.05;
edge_col_pos = [0.5835    0.6526    0.7656];
edge_col_neg = [0.9569    0.7294    0.7294];
for i = 1:(size(rho,1)-1)
    for j = (i+1):size(rho,2)
        if p(i,j) < alpha
            lw = abs(rho(i,j)) .* edge_scale;
            lc = edge_col_pos;
            if rho(i,j) < 0; lc = edge_col_neg; end
            plot(node_locs([i j],1),node_locs([i j],2),'-','LineWidth',lw,'Color',lc)
        end
    end
end

% Nodes
scatter(node_locs(isMotor,1),node_locs(isMotor,2),300,min(node_cols(isMotor,:).*1.3,1),'filled','Marker','s')
scatter(node_locs(~isMotor,1),node_locs(~isMotor,2),200,node_cols(~isMotor,:),'filled','Marker','o')
scatter(node_locs(isMotor,1),node_locs(isMotor,2),300,node_cols(isMotor,:).*0.7,'Marker','s','LineWidth',2)
scatter(node_locs(~isMotor,1),node_locs(~isMotor,2),200,node_cols(~isMotor,:).*0.7,'Marker','o','LineWidth',2)
% scatter(node_locs(:,1),node_locs(:,2),80,node_cols,'filled')
% scatter(node_locs(:,1),node_locs(:,2),80,node_cols.*0.7,'LineWidth',1)

% Node numbers
for i = 1:numel(labels)
    text(node_locs(i,1),node_locs(i,2),num2str(i),'Color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontName','Arial')
end

% Number legend
base_y = 550;
assay = 'EPM';
idx = find(contains(labels,assay));
text(0,base_y,af(@(i)sprintf('%d. %s',i,labels{i}),idx),'units','pixels','VerticalAlignment','top','FontSize',8,'interpreter','none','Color',cm(1,:),'FontWeight','bold')

assay = 'SC';
idx = find(contains(labels,assay));
text(0,base_y-((idx(1)-1)*15),af(@(i)sprintf('%d. %s',i,labels{i}),idx),'units','pixels','VerticalAlignment','top','FontSize',8,'interpreter','none','Color',cm(2,:),'FontWeight','bold')

assay = 'YM';
idx = find(contains(labels,assay));
text(0,base_y-((idx(1)-1)*15),af(@(i)sprintf('%d. %s',i,labels{i}),idx),'units','pixels','VerticalAlignment','top','FontSize',8,'interpreter','none','Color',cm(3,:),'FontWeight','bold')

% Legend
h_motor = scatter(NaN,NaN,300,[0 0 0]+0.3,'Marker','s','LineWidth',2);
h_nonmotor = scatter(NaN,NaN,200,[0 0 0]+0.3,'Marker','o','LineWidth',2);
h_sep = plot(NaN,NaN,'.','Color','w');
edge_cols = {edge_col_neg, edge_col_pos};
edge_leg_rhos = [-1 -0.5 -0.25 0.25 0.5 1];
h_edges = arrayfun(@(x)plot(NaN,NaN,'-','LineWidth',abs(x)*edge_scale,'Color',edge_cols{(x<0)+1}),edge_leg_rhos);
legend([h_motor,h_nonmotor,h_sep,h_edges],[{'Motor','Nonmotor',''}, af(@(x)sprintf('\\rho = %.2f',x), edge_leg_rhos)])

% Formatting
axis equal
noax
painters
end
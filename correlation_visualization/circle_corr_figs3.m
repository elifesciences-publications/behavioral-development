% Goal here is to make circle corr figs with sig diffs
cc
%% Setup
WT_corrPath = 'data/WT.xlsx';
corrPaths = {'data/Adult_Crus I.xlsx'
'data/Adult_Crus II.xlsx'
'data/Adult_Lobule VI.xlsx'
'data/Adult_Lobule VII.xlsx'
'data/Juven_Crus I.xlsx'
'data/Juven_Crus II.xlsx'
'data/Juven_Lobule VI.xlsx'
'data/Juven_Lobule VII.xlsx'};

alpha = 0.05;

%%
motor_labels = {
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
isMotor = cell2mat(motor_labels(:,2));

% motor_labels = isMotor([isMotor{:,2}],1);
% nonmotor_labels = isMotor(~[isMotor{:,2}],1);

%%
labels = motor_labels(:,1);
[grps,names] = findgroups(extractBefore(labels,'_'));
N = arrayfun(@(i)sum(grps==i),unique(grps));

theta = af(@(i)linspace(0,2*pi,N(i)+1),unique(grps));
theta = af(@(x)theta{x} + ((pi/4) * (x-0)),1:numel(theta));
xd = af(@(i)cos(theta{i}(1:end-1)),unique(grps));
yd = af(@(i)sin(theta{i}(1:end-1)),unique(grps));

xd{1} = xd{1} - 3; yd{1} = yd{1} + 0;
xd{2} = xd{2} + 3; yd{2} = yd{2} + 0;
xd{3} = xd{3} + 0; yd{3} = yd{3} + 3;

xd = cellcat(xd,2)'; yd = cellcat(yd,2)';
xyd = [xd yd];

%%
WT_Z = readtable(WT_corrPath,'Sheet','Z','ReadRowNames',true);
WT_rho = readtable(WT_corrPath,'Sheet','rho','ReadRowNames',true);
WT_p = readtable(WT_corrPath,'Sheet','p','ReadRowNames',true);
WT_isSig = WT_p{:,:} < 0.05;

%%
for c = 1:numel(corrPaths)
corrPath = corrPaths{c};
Z = readtable(corrPath,'Sheet','Z','ReadRowNames',true);
rho = readtable(corrPath,'Sheet','rho','ReadRowNames',true);
p = readtable(corrPath,'Sheet','p','ReadRowNames',true);
isSig = p{:,:} < alpha;

%%
alpha = 0.05;

rho_scale = 10;

cm = viridis(2);
motor_col = cm(1,:);
nonmotor_col = cm(2,:);

edge_col = [0.7294    0.8157    0.9569] .* 0.8;
neg_col = [0.9569    0.7294    0.7294];

% wt_edge_col(4) = 0.5;
% edge_col(4) = 0.5;
% xor_edge_col(4) = 0.9;
% neg_col(4) = 0.9;

figure,figclosekey
figsize(1000,750)
hold on

for i = 1:numel(labels)
    for j = (i+1):numel(labels)
        WT_rho_ij = WT_rho{i,j};
        WT_p_ij = WT_p{i,j};
        WT_isSig_ij = WT_p_ij < alpha;
        
        rho_ij = rho{i,j};
        p_ij = p{i,j};
        isSig_ij = (p_ij < alpha & sign(WT_rho_ij) ~= sign(rho_ij)) | (p_ij < alpha & ~WT_isSig_ij);
%         isSig_ij = p_ij < alpha; % all_sig

        if isSig_ij
            col = edge_col;
            if rho_ij < 0; col = neg_col; end
            plotpts(xyd([i,j],:),'-','color',col,'linewidth',rho_scale .^ abs(rho_ij))
        end
    end
end

hnodes = {};
hnodes{1} = plotpts(xyd(isMotor,:),'o','markersize',20,'markerfacecolor',motor_col,'markeredgecolor',motor_col .* 0.7);
hnodes{2} = plotpts(xyd(~isMotor,:),'o','markersize',20,'markerfacecolor',nonmotor_col,'markeredgecolor',nonmotor_col .* 0.7);
hnodes = cellcat(hnodes);

legend(hnodes,{'Motor','Non-motor'},'fontsize',20)
legend('boxoff')

for i = 1:numel(labels)
    text(xyd(i,1),xyd(i,2),num2str(i),'Color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontName','Arial')
end
text(0,1,af(@(i)sprintf('%d. %s',i,labels{i}),1:numel(labels)),'units','normalized','VerticalAlignment','top','FontSize',8,'interpreter','none')

title(replace(get_filename(corrPath,true),'_','/'),'FontSize',18)
noax

export_fig(['figs/v2_alpha=0.05/' get_filename(corrPath,true)],'-depsc','-painters','-eps','-png')
% export_fig(['figs/v2_alpha=0.05/all_sig' get_filename(corrPath,true)],'-depsc','-painters','-eps','-png')
close(gcf)

end

%%
alpha = 0.05;

rho_scale = 10;

cm = viridis(2);
motor_col = cm(1,:);
nonmotor_col = cm(2,:);

edge_col = [0.7294    0.8157    0.9569] .* 0.8;
neg_col = [0.9569    0.7294    0.7294];

% edge_col(4) = 0.5;
% xor_edge_col(4) = 0.9;

figure,figclosekey
figsize(1000,750)
hold on
for i = 1:numel(labels)
    for j = (i+1):numel(labels)
%         WT_rho_ij = WT_rho{i,j};
%         WT_p_ij = WT_p{i,j};
%         WT_isSig_ij = WT_p_ij < alpha;
        
        rho_ij = WT_rho{i,j};
        p_ij = WT_p{i,j};
        isSig_ij = p_ij < alpha;

        if isSig_ij
            col = edge_col;
            if rho_ij < 0; col = neg_col; end
            plotpts(xyd([i,j],:),'-','color',col,'linewidth',rho_scale .^ abs(rho_ij))
        end
    end
end

hnodes = {};
hnodes{1} = plotpts(xyd(isMotor,:),'o','markersize',20,'markerfacecolor',motor_col,'markeredgecolor',motor_col .* 0.7);
hnodes{2} = plotpts(xyd(~isMotor,:),'o','markersize',20,'markerfacecolor',nonmotor_col,'markeredgecolor',nonmotor_col .* 0.7);
hnodes = cellcat(hnodes);

legend(hnodes,{'Motor','Non-motor'},'fontsize',20)
legend('boxoff')

for i = 1:numel(labels)
    text(xyd(i,1),xyd(i,2),num2str(i),'Color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontName','Arial')
end
text(0,1,af(@(i)sprintf('%d. %s',i,labels{i}),1:numel(labels)),'units','normalized','VerticalAlignment','top','FontSize',8,'interpreter','none')

title('WT','FontSize',18)
noax

export_fig('figs/v2_alpha=0.05/WT','-depsc','-painters','-eps','-png')
close(gcf)

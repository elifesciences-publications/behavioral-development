%% Setup
% Let's start clean :)
clear, clc, close all

% Add dependencies to the path
addpath(genpath('deps'))

% Load data
load('data/2017-09-04-master_dataset.mat')

%% Adult vs Juven
alpha = 0.000001;

figDir = 'figs/adult_vs_juven_2d';
for r = 1:numel(regions)
    region = regions{r};
    
    % Setup saving
    saveDir = ff(figDir, region);
    mkdirto(saveDir)
    
    % Pull out region data
    X_region = X(behavior.Region(X.MouseID) == string(region),:);
    
    % Pull out age groups
    X1 = X_region(behavior.Age(X_region.MouseID) == "Adult",:);
    X2 = X_region(behavior.Age(X_region.MouseID) == "Juven",:);

    % Initialize container to the number of metrics^2
    D = size(X_region,2);
    p = NaN(D);
    
    % Iterate through every pair of metrics
    for i = 1:D
        for j = 1:D
            % Pull out data for pair of metrics i and j
            A = rmmissing(X1{:,[i j]});
            B = rmmissing(X2{:,[i j]});
            
            % Run 2d K-S test
            test = DepTest2(A,B,'test','ks');
            p(i,j) = test.pval;
            
            if p(i,j) >= alpha; continue; end
            
            % Plotting
            labels = {sprintf('Acute %s (n = %d)', region, size(A,1)) 
                sprintf('Developmental %s (n = %d)', region, size(B,1))};
            if i == j
                % Plot boxplot
                fig = figure;
                [X_i, G] = cellcat({A(:,1);B(:,1)});
                boxplot(X_i,G,'Orientation','horizontal')
                h = plotSpread(X_i,'distributionIdx',G,'distributionColors',{'r','b'},'xyOri','flipped');
                set(h{1},'markersize',25)
                boxplotformat
                xlabel(metrics{i},'Interpreter','none')
                yticks([1 2])
                yticklabels(labels)
                graygrid
                fontsize(16)
                painters
                
                figsize(1100,600)
            else
                % Plot scatter
                fig = figure();
                h = [];
                h(1) = plotpts(A,'r.'); hold on
                h(2) = plotpts(B,'b.');
                set(h, 'MarkerSize', 20)
                graygrid
                xlabel(metrics{i},'Interpreter','none','FontSize',16)
                ylabel(metrics{j},'Interpreter','none','FontSize',16)
                legend(labels,'Color','w','Location','northeastoutside')
                fontsize(16)
                painters
                
                figsize(1100,600)
                axis equal
            end
            titlef('p = %f (K-S; null: same distribution)', p(i,j))
            
            % Save
            savePath = ff(saveDir, sprintf('%s-%s', metrics{i}, metrics{j}));
            export_fig(savePath, '-png', '-eps')
            close(gcf)
        end
    end
    
    figure, figclosekey
    histogram(p)
    xlabel('p')
    title(region)
    
    continue
    % Plot region summary
    figure, figclosekey
    imagesc(p)
    colorbar2('p-value')
    axis tight ij
    colormap(viridis)
    xlabel('Developmental','FontSize',14), ylabel('Adult','FontSize',14)
    xticks(1:size(p,2)), xticklabels(metrics), xtickangle(90)
    yticks(1:size(p,1)), yticklabels(metrics)
    shortticks
    ax = gca;
    ax.TickLabelInterpreter = 'none';
    colormap(flipud(inferno))
    % colormap(flipud(jet))
    figsize(1100,900)
    title(region,'FontSize',18)
    savePath = ff(saveDir, '_summary');
    export_fig(savePath, '-png', '-eps')
    close(gcf)
end

%% Experimental vs control
alpha = 0.01;

figDir = 'figs/expt_vs_control_2d';
for a = 1:numel(ages)
    age = ages{a};
    
    % Pull out age data
    M_age = M(:,:,a);
    
    for r = 1:numel(regions)
        region = regions{r};

        % Setup saving
        saveDir = ff(figDir, [age ' ' region]);
        mkdirto(saveDir)

        % Pull out region data
        M_region = M_age(r,:);
        M_region = af(@(m) M_region{m}(:,metrics{m}), 1:numel(metrics));
        
        % Pull out groups
%         X1 = X_region(behavior.Age(X_region.MouseID) == "Adult",:);
%         X2 = X_region(behavior.Age(X_region.MouseID) == "Juven",:);

        % Initialize container to the number of metrics^2
        D = numel(metrics);
        p = NaN(D);

        % Iterate through every pair of metrics
        for i = 1:D
            for j = 1:D
                % Pull out data for mice in both
                mice = intersect(M_region{i}.MouseID, M_region{j}.MouseID);
                M_ij = [M_region{i}{mice,:}, M_region{j}{mice,:}];
                
                % Separate by condition
                A = rmmissing(M_ij(behavior.Condition(mice) == "Control",:));
                B = rmmissing(M_ij(behavior.Condition(mice) == "Experimental",:));

                % Run 2d K-S test
                test = DepTest2(A,B,'test','ks');
                p(i,j) = test.pval;

                if p(i,j) >= alpha; continue; end

                % Plotting
                labels = {sprintf('%s %s Control (n = %d)', age, region, size(A,1)) 
                    sprintf('%s %s Experimental (n = %d)', age, region, size(B,1))};
                if i == j
                    % Plot boxplot
                    fig = figure;
                    [X_i, G] = cellcat({A(:,1);B(:,1)});
                    boxplot(X_i,G,'Orientation','horizontal')
                    h = plotSpread(X_i,'distributionIdx',G,'distributionColors',{'r','b'},'xyOri','flipped');
                    set(h{1},'markersize',25)
                    boxplotformat
                    xlabel(metrics{i},'Interpreter','none')
                    yticks([1 2])
                    yticklabels(labels)
                    graygrid
                    fontsize(16)
                    painters

                    figsize(1100,600)
                else
                    % Plot scatter
                    fig = figure();
                    h = [];
                    h(1) = plotpts(A,'r.'); hold on
                    h(2) = plotpts(B,'b.');
                    set(h, 'MarkerSize', 20)
                    graygrid
                    xlabel(metrics{i},'Interpreter','none','FontSize',16)
                    ylabel(metrics{j},'Interpreter','none','FontSize',16)
                    legend(labels,'Color','w','Location','northeastoutside')
                    fontsize(16)
                    painters

                    figsize(1100,600)
                    axis equal
                end
                titlef('p = %f (K-S; null: same distribution)', p(i,j))

                % Save
                savePath = ff(saveDir, sprintf('%s-%s', metrics{i}, metrics{j}));
                export_fig(savePath, '-png', '-eps')
                close(gcf)
            end
        end

        % Plot region summary
        figure, figclosekey
        imagesc(p)
        colorbar2('p-value')
        axis tight ij
        colormap(viridis)
%         xlabel('Control','FontSize',14), ylabel('Experimental','FontSize',14)
        xticks(1:size(p,2)), xticklabels(metrics), xtickangle(90)
        yticks(1:size(p,1)), yticklabels(metrics)
        shortticks
        ax = gca;
        ax.TickLabelInterpreter = 'none';
        colormap(flipud(inferno))
        % colormap(flipud(jet))
        figsize(1100,900)
        title([age ' ' region],'FontSize',18)
        savePath = ff(saveDir, '_summary');
        export_fig(savePath, '-png', '-eps')
        close(gcf)
    end
end
% This script will run an analysis on the behavior data treating each assay
% individually for every age and region, computing the PCA with respect to
% the specified controls and then performing a two-sample
% Kolmogorov-Smirnov test to see which metrics were likely drawn from the 
% same distribution.
%
% Originally in: ../2016-04-12 - behavior pca - the separability awakens

%% Setup
% Let's start clean :)
clear, clc, close all

% Add dependencies to the path
addpath(genpath('deps'))

%% General parameters
% These will apply to all conditions that will be compared

% Path to data spreadsheet
% dataPath = 'BehavioralDataMatrix.April14.xlsx';
dataPath = 'data/BehavioralDataMatrix.April14_clean_8_17_2017.xlsx';

% Path to where data will be saved
% savePath = 'data/ks_analysis.mat';
savePath = ['data', filesep, isodate(), '-ks_analysis.mat'];
% savePath = ''; % nothing will be saved

% Experimental regions
exptRegions = {
    'Crus I'
    'Crus II'
    'Lobule VI'
    'Lobule VII'
    };

% Experimental ages
exptAges = {
    'Adult'
    'Juven'
%     'Infant'
    };

% Completeness of behavior table
defaults.completeness = 'metric'; % keep any mouse that has data for the specified metrics
% defaults.completeness = 'all';    % exclude all mice missing ANY data

% Mice to exclude from the entire analysis
defaults.excludeMice = {
    % Weird outliers
    'GR3_2'
    'GR3W7_4'
    % Excitatory DREADDs
    'GR15_1'
    'GR15_2'
    'GR15_3'
    'GR15_13'
    'GR15_14'
    'GR15_5'
    'GR15_6'
    'GR15_8'
    'GR15_10'
    'GR15_11'
    'GR15_12'
    'GR15_4'
    'GR15_7'
    'GR15_9'
    'GR15_15'
    % Bad early mice
    'GRBL6T_1'
    'GRBL6T_8'
    'GRBL6S_3'
    'GRBL6S_4'
    'GRBL6S_6'
    };

%% Per-assay configuration
% Y-maze
ym.name = 'Ymaze';
ym.abbrev = 'YM';
ym.desc = 'Y-maze';
ym.params = defaults;
ym.params.keepMetrics = {
    'YM_AcqInitialLR'
    'YM_AcqSecondaryLR'
    'YM_AcqAbility'
    'YM_EarlyRevInitialLR'
    'YM_EarlyRevSecondaryLR'
    'YM_EarlyRevAbility'
    'YM_LateRevInitialLR'
    'YM_LateRevAbility'
    'YM_Distance'
    'YM_Velocity'
};
ym.params.cnoAges = exptAges; % include CNO controls from all ages
ym.params.includeWT = true;

% Grooming
gr.name = 'Grooming';
gr.abbrev = 'GR';
gr.desc = 'Grooming';
gr.params = defaults;
gr.params.keepMetrics = {'GroomingRatio'};
gr.params.cnoAges = {};
gr.params.includeWT = false;

% Social chamber
sc.name = 'SocialChamber';
sc.abbrev = 'SC';
sc.desc = 'Social Chamber';
sc.params = defaults;
sc.params.keepMetrics = {
    'SC_NoveltySeeking'
    'SC_SocialPreference'
    'SC_BSDistance'
    'SC_TestDistance'
    };
sc.params.cnoAges = exptAges; % include CNO controls from all ages
sc.params.includeWT = true;

% EPM
epm.name = 'EPM';
epm.abbrev = 'EPM';
epm.desc = 'EPM'; % 'Elevated Plus Maze'
epm.params = defaults;
epm.params.keepMetrics = {
    'EPM_OpenArmPref'
    'EPM_AntiHesitation'
    'EPM_ExplorationTime'
    'EPM_ExplorationEntrances'
    'EPM_Distance'
    'EPM_Velocity'
    };
epm.params.cnoAges = {}; % include only age-matched CNO controls
epm.params.includeWT = false;

%% Process
warning('off', 'stats:pca:ColRankDefX')
assays = {ym, gr, sc, epm};
kstable = table();
for i = 1:numel(assays)
    assay = assays{i};
    for j = 1:numel(exptAges)
        exptAge = exptAges{j};
        for k = 1:numel(exptRegions)
            exptRegion = exptRegions{k};
            printf('Processing: *%s* -> *%s* -> *%s*\n\t #nonl', assay.name, exptAge, exptRegion)

            % Load and filter behavior data table
            B = loadBehavior(dataPath, assay.params, 'keepAges', {exptAge}, ...
                    'keepRegions', {exptRegion});

            % PCA
            B.pcs = pcaproj(B.metrics, B.isCtrl, ~B.isCtrl);
            
            % Assay naming consistency
            B.assay_name = assay.name;
            B.assay_abbrev = assay.abbrev;
            noAbbrev = ~startsWith(B.metricLabels,[assay.abbrev, '_']);
            B.metricLabels(noAbbrev) = strcat(assay.abbrev, '_', B.metricLabels(noAbbrev));
            
            % Define some tables for convenience
            %   Raw metrics
            B.raw = array2table(B.metrics, 'Var', B.metricLabels);
            B.raw.Properties.VariableDescriptions = B.descriptions;
            B.raw.Properties.VariableUnits = B.units;
            %   Normalized metrics
            B.X = array2table(B.pcs.X, 'Var', B.metricLabels);
            B.X.Properties.VariableDescriptions = B.descriptions;
            B.X.Properties.VariableUnits = repmat({'Control Standardized'}, size(B.metricLabels));
            %   PCA
            B.pca = array2table(B.pcs.score, 'Var', B.pcs.labels);
            B.pca.Properties.VariableNames = strcat(assay.abbrev, '_', B.pca.Properties.VariableNames);
            B.pca.Properties.VariableDescriptions = strcat(assay.desc, {' '}, B.pcs.labels);
            B.pca.Properties.VariableUnits = repmat({'PC Score'}, size(B.pcs.labels));
            %   Normalized metrics and PCA
            B.all = [B.X, B.pca];
            % Make table indexable by mouse name like .behavior
            B.all.Properties.RowNames = B.mice;
            %   Split by condition
            B.expt = B.all(~B.isCtrl,:);
            B.ctrl = B.all(B.isCtrl,:);
            %   Convenient access to text and captions
            B.allLabels = B.all.Properties.VariableNames;
            B.allDescriptions = B.all.Properties.VariableDescriptions;
            B.allUnits = B.all.Properties.VariableUnits;
            B.rawUnit = cell2struct(B.units(:), B.metricLabels(:));
            B.desc = cell2struct(B.allDescriptions(:), B.allLabels(:));
            B.unit = cell2struct(B.allUnits(:), B.allLabels(:));
            
            % Compute KS table
            for l = 1:numel(B.allLabels)
                label = B.allLabels{l};
                [~, B.ks.(label)] = kstest2(B.expt{:,label}, B.ctrl{:,label});
            end
            
            % More metadata
            regionVarName = strrep(exptRegion, ' ', '');
            B.region_name = regionVarName;
            B.region_desc = exptRegion;
            B.age_name = exptAge;
            
            % Save processing results
            assay.(regionVarName).(exptAge) = B;
            
            % Append observations to KS table
            obs = table();
            obs{:,'VarName'} = vert(fieldnames(B.ks));
            obs{:,'p'} = vert(struct2array(B.ks));
            obs{:,'Assay'} = {assay.name};
            obs{:,'Age'} = {exptAge};
            obs{:,'Region'} = {exptRegion};
            obs{:,'NumExpt'} = sum(~B.isCtrl);
            obs{:,'NumCtrl'} = sum(B.isCtrl);
            obs{:,'DeltaMu'} = vert(varfun(@mean, B.expt, 'Out', 'uni') - varfun(@mean, B.ctrl, 'Out', 'uni'));
            obs{:,'Data'} = {B}; % extreme redundancy for convenience
            kstable = [kstable; obs];
        end
    end
    assays{i} = assay;
end
[ym, gr, sc, epm] = assays{:};
kstable = kstable(:, {'Region', 'Age', 'Assay', 'VarName', 'p', 'NumExpt', 'NumCtrl', 'DeltaMu', 'Data'});
warning('on', 'stats:pca:ColRankDefX')

% Save
if ~isempty(savePath)
    savePath = get_new_filename(savePath);
    timestamp = datestr(now);
    save(savePath, 'ym', 'gr', 'sc', 'epm', 'kstable', 'exptRegions', 'exptAges', 'timestamp')
    printf('Saved results to: *%s*', savePath)
end

printf
%% Display significant
alpha = 0.05;
minDelta = 0.1;

isSignificant = kstable.p < alpha & abs(kstable.DeltaMu) > minDelta;
significant = sortrows(kstable(isSignificant, :), {'Region', 'Age', 'Assay', 'p'});

printf('*Significant observations* (alpha = %g, minDeltaMu = %g):', alpha, minDelta)
disp(significant)

return

%% Visualization
figure, plotRaw(epm.LobuleVII.Juven, 'EPM_OpenArmPref')
figure, plotMetric(sc.CrusI.Juven, 'SC_SocialPreference')
figure, plotPC(sc.CrusI.Juven, 3)
figure, plotPCA(sc.CrusI.Juven)

% This script takes the KS analysis data and compiles it into a more useful
% format for further analysis.
%
% This is similar to the masterlist but with more redundancies.
%
% Check the bottom of this file for descriptions of all of the data that is
% saved here.

%% Setup
% Let's start clean :)
clear, clc, close all

% Add dependencies to the path
addpath(genpath('deps'))

% Raw behavior metrics
rawDataPath = 'data/BehavioralDataMatrix.April14_clean_8_17_2017.xlsx';

% Analysis data generated by ks_analysis.m
ksAnalysisPath = 'data/2017-08-31-ks_analysis.mat';

% Save path for merged dataset
masterDatasetPath = ['data/' isodate() '-master_dataset.mat'];

%% Load raw data
% Read and parse data
opts = detectImportOptions(rawDataPath, 'NumHeaderLines',2);
behavior = readtable(rawDataPath, opts);

% Units
[~,~,raw] = xlsread(rawDataPath);
behavior.Properties.VariableDescriptions = raw(1,:);
units = raw(2,:); units(~cellfun(@(x)ischar(x),units)) = {''};
behavior.Properties.VariableUnits = units;

% Make indexable by mouse ID
behavior.Properties.RowNames = behavior.MouseID;
behavior.Properties.DimensionNames = {'MouseID','Metrics'};
clear raw units opts

%% Load KS analysis data
load(ksAnalysisPath)

%% Create merged experimental table
data = cellcat(kstable.Data);

% Find all experimental mice and metrics
X_mice = unique(cellcat(af(@(x)x.mice(~x.isCtrl), data)));
X_metrics = unique(cellcat(af(@(x)x.allLabels, data),2));

% Create experimental data table
X = array2table(NaN(numel(X_mice),numel(X_metrics)), 'RowNames', X_mice, 'VariableNames', X_metrics);

% Fill in merged experimental table
for i = 1:numel(data)
    r = data(i).expt.Properties.RowNames;
    c = data(i).expt.Properties.VariableNames;
    X(r,c) = data(i).expt;
end
X.Properties.DimensionNames = {'MouseID','Metrics'};

% Gather metadata
X_meta = behavior(X.Properties.RowNames,{'MouseID','Age','Region','Condition'});

%% Create experiment matrix: regions x metrics x age (like SfN figs)
% Find dimension names
regions = unique(X_meta.Region);
metrics = X_metrics;
ages = unique(X_meta.Age);

% Initialize
M = cell(numel(regions), numel(metrics), numel(ages));
M_pca = M;

% Fill in the matrix
for i = 1:height(kstable)
    % Find where we are in the matrix
    r = strcmp(regions, kstable.Region{i});
    m = strcmp(metrics, kstable.VarName{i});
    a = strcmp(ages, kstable.Age{i});
    
    % Pull out normalized data
    M{r,m,a} = kstable.Data{i}.all;
    M{r,m,a}.Properties.DimensionNames = {'MouseID', 'Metrics'};
    
    % Pull out PCA metadata
    M_pca{r,m,a} = kstable.Data{i}.pcs;
end

%% Save
timestamp = datestr(now);
vars = {
    'timestamp'
    
    % Paths
    'rawDataPath' % raw data spreadsheet
    'ksAnalysisPath' % from ks_analysis.m
    'masterDatasetPath' % output of this script
    
    % Parsed behavior matrix spreadsheet
    'behavior' % table with all of the raw data
    
    % Experimental mice merged table
    'X_mice' % all experimental mice
    'X_metrics' % all metric names
    'X' % all of the normalized behavior metrics for the experimental mice
    'X_meta' % metadata for X (mouseID, age, region, condition)
    
    % Experiment matrix (like SfN figs)
    'regions' % all of the anatomical regions
    'metrics' % all of the metrics (same as X_metrics)
    'ages' % all of the mouse ages (i.e., {'Adult', 'Juven'})
    'M' % cell array of (regions x metrics x ages) containing both experimental and matched control data
    'M_pca' % the PCA data (coeffs, latents, etc) that matches up with M
    };
save(masterDatasetPath, vars{:})

% This script merges the registered anatomical expression data into a
% single file for analysis.
% 
% Check the bottom of this file for descriptions of all of the data that is
% saved here.

%% Setup
% Let's start clean :)
clear, clc, close all

% Add dependencies to the path
addpath(genpath('deps'))

% Clean signal stacks folder (data in cleaned_tiffs.zip)
% stacksPath = 'data/cleaned_tiffs';
stacksPath = '\\bucket.pni.princeton.edu\wang\Talmo\data\DREADDs\TV\cleaned_tiffs';

% Cropped and annotated Allen Brain Atlas
atlasStackPath = 'data/atlas/avg_25um.coronal.tif';
atlasMeshPath = 'data/atlas/mesh.mat';
bilateral_rois = {'Crus 1','Crus 2','Paramedian lobule','Simplex lobule'};

% Save path for merged anatomical data
anatomyPath = ['data/' isodate() '-anatomy.mat'];

%% Load and merge atlas data
% Load image data
vol = readTIFFStack(atlasStackPath);

% Pull out region stacks
atlas = load(atlasMeshPath);

% Split bilateral ROIs into left/right
BW = vert({atlas.contours.BW});
rois = {atlas.contours.type};
for i = 1:numel(bilateral_rois)
    idx = strcmpi(rois,bilateral_rois{i});
    roi = atlas.contours(idx);
    BWi = roi.BW;
    
    BW_right = BWi; BW_right(:,1:end/2,:) = 0;
    BW{end+1} = BW_right;
    rois{end+1} = [roi.type ' (Right)'];
    
    BW_left = BWi; BW_left(:,end/2+1:end,:) = 0;
    BW{end+1} = BW_left;
    rois{end+1} = [roi.type ' (Left)'];
    
end
for i = 1:numel(bilateral_rois)
    idx = strcmpi(rois,bilateral_rois{i});
    BW(idx) = [];
    rois(idx) = [];
end
[rois,idx] = natsort(rois);
BW = BW(idx);
% BW = cellcat(BW,4);

% Merge
atlas = varstruct(atlasStackPath, vol, atlasMeshPath, rois, BW);

%% Load sample expression data
% Find sample data
samplePaths = searchFor('*.tif', stacksPath);
samples = get_filename(samplePaths,true); % sample names from filenames

% Load all the signal data
signal = cf(@(x)readTIFFStack(x) > 0, samplePaths);

% Merge into a single array
% signal = cellcat(signal, 4);

%% Compute ROI coverage statistics
% Initialize tables
signal_distribution = NaN(numel(samples), numel(atlas.rois));
roi_coverage = NaN(numel(samples), numel(atlas.rois));
jaccard_idx = NaN(numel(samples), numel(atlas.rois));

% Fill in table
for i = 1:numel(samples)
% for i = 1
    % Pull out and count signal voxels for sample
    X = signal{i};
    nX = sum(X(:));
    
    for j = 1:numel(atlas.rois)
        % Pull out and count voxels for ROI
        Y = atlas.BW{j};
        nY = sum(Y(:));
        
        % Compute intersection and union
        nXxY = sum(X(:) & Y(:));
        nXuY = sum(X(:) | Y(:));
        
        % Compute metrics
        signal_distribution(i,j) = nXxY / nX;
        roi_coverage(i,j) = nXxY / nY;
        jaccard_idx(i,j) = nXxY / nXuY;
    end
end

% Sanitize ROI names
roi_vars = replace(titlecase(atlas.rois), {' ','(',')','/'}, '');

% Convert to tables
signal_distribution = array2table(signal_distribution, 'RowNames', samples, 'VariableNames', roi_vars);
roi_coverage = array2table(roi_coverage, 'RowNames', samples, 'VariableNames', roi_vars);
jaccard_idx = array2table(jaccard_idx, 'RowNames', samples, 'VariableNames', roi_vars);

% Add ROI labels
signal_distribution.Properties.VariableDescriptions = atlas.rois;
roi_coverage.Properties.VariableDescriptions = atlas.rois;
jaccard_idx.Properties.VariableDescriptions = atlas.rois;

%% Save
timestamp = datestr(now);
vars = {
    'timestamp'
    
    % Paths
    'stacksPath' % path to folder containing clean sample signal TIFFs (cleaned_tiffs.zip)
    'atlasStackPath' % raw Atlas TIFF intensity stack
    'atlasMeshPath' % preprocessed meshed Atlas
    'bilateral_rois' % names of ROIs that are bilateral
    'anatomyPath' % output of this script
    
    % Atlas
    'atlas' % struct with Allen Brain Atlas data:
            %   .vol: raw intensity stack (coronal)
            %   .rois: names of each ROI
            %   .BW: cell array with each ROI as a binary stack of the same
            %        size as vol
    
    % Samples (registered to Allen Brain Atlas)
    'samples' % names for mice with expression data
    'signal' % cell array with a binary stack for each sample
    'samplePaths' % paths to each sample signal stack
    
    % Coverage statistics
    'signal_distribution' % fraction of sample signal in each ROI
    'roi_coverage' % fraction of each ROI covered by signal
    'jaccard_idx' % Jaccard index (intersection / union)
    };
save(anatomyPath, vars{:})

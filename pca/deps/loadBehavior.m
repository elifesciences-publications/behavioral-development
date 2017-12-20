function X = loadBehavior(dataPath, varargin)

% Set parameter defaults
defaults = struct();
defaults.xlsRange = ''; % autosize
defaults.keepMetrics = {};
defaults.keepAges = {};
defaults.keepRegions = {};
defaults.excludeMice = {};
defaults.cnoAges = {}; % empty = match expt ages
defaults.includeWT = true; % true = include WT controls
defaults.completeness = 'all'; % 'all' or 'metric'
defaults.metafields = {'Age', 'Region', 'Condition'};
defaults.IDfield = 'MouseID';
defaults.verbosity = 0;

% Spreadsheet format
descRow = 1;
unitsRow = 2;
varsRow = 3;
dataStartRow = 4;

% Parse input parameters
params = parse_params(varargin, defaults);

% Load spreadsheet
[~, ~, raw] = xlsread(dataPath, params.xlsRange);

% Parse out variable names
varNames = cf(@strtrim, raw(varsRow, :));
varNames = cf(@(x) strrep(x, ' ', '_'), varNames);
data = raw(dataStartRow:end, :);

% Parse out metafields
isMetafield = ismember(varNames, params.metafields);
metafields = varNames(isMetafield);

% Parse out ID field (mouse ID)
isIDfield = strcmp(varNames, params.IDfield);
metricLabels = varNames(~isIDfield & ~isMetafield);

% Parse out descriptions
descriptions = raw(descRow, :);
descriptions(cellfun(@(x) all(isnan(x)), descriptions)) = {''};

% Parse out units
units = raw(unitsRow, :);
units(cellfun(@(x) all(isnan(x)), units)) = {''};

% Enable indexing with mouse names as well
mice = data(:, isIDfield);

% Construct table
behavior = cell2table(data, 'Var', varNames, 'Row', mice);
behavior.Properties.VariableDescriptions = descriptions;
behavior.Properties.VariableUnits = units;

% Start off including everyone
idx = true(height(behavior), 1);

% Keep only mice of the specified ages
if ~isempty(params.keepAges)
    isKeptAge = ismember(behavior.Age, params.keepAges);
    idx = idx & isKeptAge;
end

% Keep only mice of the specified regions
if ~isempty(params.keepRegions)
    isKeptRegion = ismember(behavior.Region, params.keepRegions);
    idx = idx & isKeptRegion;
end

% Filter CNO controls by age
isCtrl = strcmpi(behavior.Condition, 'Control');
isCNOCtrl = strcmp(behavior.Region, 'CNO Ctl');
if isempty(params.cnoAges)
    params.cnoAges = unique(behavior.Age(idx));
end
isCNOAge = ismember(behavior.Age, params.cnoAges) & isCNOCtrl;
idx = (idx & ~isCtrl) | isCNOAge;

% Include wildtype controls
if params.includeWT
    isWT = strcmp(behavior.Region, 'WT Control');
    idx = idx | isWT;
end

% Exclude specific mice
if ~isempty(params.excludeMice)
    excludedMice = ismember(behavior.MouseID, params.excludeMice);
    idx = idx & ~excludedMice;
end

% Exclude any mice that are missing data for ANY metric
if strcmp(params.completeness, 'all')
    missingData = any(ismissing(behavior(idx, metricLabels)), 2);
    idx = idx & ~missingData;
end

% Keep only specified metrics
if ~isempty(params.keepMetrics)
    assert(all(ismember(params.keepMetrics, metricLabels)), ...
        'Specified metrics to keep that are not in the data.')
    metricLabels = params.keepMetrics(:)';
end

% Exclude any mice that are missing data with the specified metrics
if strcmpi(params.completeness, 'metric')
    missingMetricData = any(ismissing(behavior(:, metricLabels)), 2);
    idx = idx & ~missingMetricData;
end

% Crop data now that we're done filtering
behavior = behavior(idx, [params.IDfield, metafields, metricLabels]);
isCtrl = strcmpi(behavior.Condition, 'Control');
isWT = strcmp(behavior.Region, 'WT Control');
printf('Loaded data for *%3d/%d* mice: *%3d* expt, *%3d* CNO ctrl, *%3d* WT ctrl', ...
    height(behavior), numel(idx), sum(~isCtrl), sum(isCtrl & ~isWT), sum(isWT));

% Return a structure with everything for portability
X = save2struct('behavior', 'metricLabels', 'metafields', 'params', 'dataPath', 'isCtrl');

% Convenience
X.mice = X.behavior.MouseID;
% X.isCtrl = strcmp(X.behavior.Condition, 'Control');
X.metrics = table2array(X.behavior(:, X.metricLabels));
X.isMetric = ismember(X.behavior.Properties.VariableNames, X.metricLabels);
X.descriptions = X.behavior.Properties.VariableDescriptions(X.isMetric);
X.units = X.behavior.Properties.VariableUnits(X.isMetric);

% Sanity checking
assert(~any(any(ismissing(X.behavior))))
assert(~any(isnan(X.metrics(:))))
assert(isequal(X.behavior{:,params.IDfield}, X.behavior.Properties.RowNames))
assert(isequal(X.behavior{:,params.IDfield}, X.mice))
assert(all(strcmp(X.behavior.Condition(X.isCtrl), 'Control')))
end

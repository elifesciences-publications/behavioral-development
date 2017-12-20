function plotRaw(B, varName, plotPoints)
%PLOTRAW Plots the distribution of a raw metric.
% Usage:
%   plotRaw(B, varName)
%
% Example:
%   plotRaw(epm.LobuleVII.Juven, 'EPM_OpenArmPref')
%   plotRaw(epm.LobuleVII.Juven, 1) % B.metricLabels{1} == 'EPM_OpenArmPref'
%
% See also: plotMetric, plotPC

if nargin < 3; plotPoints = false; end

if isnumeric(varName)
    varName = B.metricLabels{varName};
end

X = B.raw.(varName);
units = B.rawUnit.(varName);
description = B.desc.(varName);
G = ~B.isCtrl;
XTicks = {
    sprintf('\\bfControl\\rm (n = %d)\\newline\\fontsize{9}\\it%s', sum(G == 0), ...
        strjoin(unique(B.behavior.Region(B.isCtrl)), ', '))
    sprintf('\\bfExperimental\\rm (n = %d)\\newline\\fontsize{9}\\it%s \\rightarrow %s', ...
        sum(G == 1), B.params.keepAges{1}, B.params.keepRegions{1})
    };
p = B.ks.(varName);

boxplot(X, G)
set(gca, 'TickLabelInterpreter', 'tex', 'XTickLabel', XTicks)
ylabel(units)
title({titlecase(description); ['\rm\fontsize{10} Kolmogorov-Smirnov: \itp\rm = ' formatnum(p)]})
boxplotformat

if plotPoints
    hold on
    plotSpread({X(~G), X(G)}, 'distributionMarkers', {'o', 'o'}, ...
        'distributionColors', {'b', 'r'}, 'spreadFcn', {'xp'}, ...
        'spreadWidth', 1.5)
end

resizeFig([380, 470])
end


function plotMetric(B, varName)
%PLOTMETRIC Plots the distribution of a normalized metric or PC.
% Usage:
%   plotMetric(B, varName)
%
% Example:
%   plotMetric(epm.LobuleVII.Juven, 'EPM_OpenArmPref')
%   plotMetric(epm.LobuleVII.Juven, 1) % B.allLabels{1} == 'EPM_OpenArmPref'
%   plotMetric(epm.LobuleVII.Juven, 'PC3')
%
% See also: plotRaw, plotPC

if isnumeric(varName)
    varName = B.allLabels{varName};
end

X = B.all.(varName);
units = B.unit.(varName);
description = B.desc.(varName);
G = ~B.isCtrl;
XTicks = {
    sprintf('\\bfControl\\rm (n = %d)\\newline\\fontsize{9}\\it%s', sum(G == 0), ...
        strjoin(unique(B.behavior.Region(B.isCtrl)), ', '))
    sprintf('\\bfExperimental\\rm (n = %d)\\newline\\fontsize{9}\\it%s \\rightarrow %s', ...
        sum(G == 1), B.params.keepAges{1}, B.params.keepRegions{1})
    };
p = B.ks.(varName);

hs = ishold;
boxplot(X, G), hold on
xl = xlim;
plot(xl, [0 0], '-', 'Color', [0.8 0.8 0.8])
set(gca, 'TickLabelInterpreter', 'tex', 'XTickLabel', XTicks)
ylabel(units)
title({titlecase(description); ['\rm\fontsize{10} Kolmogorov-Smirnov: \itp\rm = ' formatnum(p)]})
boxplotformat
if ~hs; hold off; end

resizeFig([380, 470])
end


function plotPC(B, pc)
%PLOTPC Plots the distribution of a principal component.
% Usage:
%   plotPC(B, pc)
%
% Example:
%   plotPC(epm.LobuleVII.Juven, 'PC3')
%   plotPC(sc.CrusI.Juven, 3)
%
% See also: plotRaw, plotMetric, plotPCA

if ischar(pc)
    varName = pc;
    if isfield(B, 'assay_abbrev')
        pc = strrep(varName, [B.assay_abbrev '_'], '');
        pc = str2double(strrep(pc, 'PC', ''));
    else
        pc = str2double(strrep(varName, 'PC', ''));
    end
elseif isnumeric(pc)
    varName = ['PC' num2str(pc)];
    if isfield(B, 'assay_abbrev')
        varName = [B.assay_abbrev '_' varName];
    end
end

X = B.pca.(varName);
units = B.unit.(varName);
description = titlecase(B.desc.(varName));
G = ~B.isCtrl;
GLabels = {
    sprintf('\\bfControl\\rm (n = %d)', sum(G == 0))
    sprintf('\\bfExperimental\\rm (n = %d)', sum(G == 1))
    };
ctrlStr = sprintf('\\fontsize{9}\\it%s', strjoin(unique(B.behavior.Region(B.isCtrl)), ', '));
exptStr = sprintf('\\fontsize{9}\\it%s \\rightarrow %s', B.params.keepAges{1}, B.params.keepRegions{1});
p = B.ks.(varName);

coeff = B.pcs.coeff(:, pc);
latent = [B.pcs.pcaLatent(pc) B.pcs.projLatent(pc)];
explained = [B.pcs.pcaExplained(pc) B.pcs.projExplained(pc)];

gap = 0.1; % vert, horz
marg_h = 0.1; % bottom, top
marg_w = 0.1; % left, right

figPos = get(gcf, 'Pos');
set(gcf, 'Pos', [figPos(1:2), 785, 340])

% Plot boxplots
subtightplot(2,2,1,gap,marg_h,marg_w)
boxplot(X, G, 'Orientation', 'horizontal', 'Label', GLabels), hold on
axis ij
yl = ylim;
plot([0 0], yl, '-', 'Color', [0.8 0.8 0.8])
set(gca, 'TickLabelInterpreter', 'tex')
xlabel(titlecase(units))
title({titlecase(description); ['\rm\fontsize{10} Kolmogorov-Smirnov: \itp\rm = ' formatnum(p)]})
tl = gca;
boxplotformat

text(-0.025, 0.525, ctrlStr, 'HorizontalAlignment', 'right', 'Units', 'normalized')
text(-0.025, 0.05, exptStr, 'HorizontalAlignment', 'right', 'Units', 'normalized')


% Plot coeffs
subtightplot(2,2,[2 4],gap,marg_h,[0.25 0.01])
plotCoeffs(coeff, B.metricLabels)
set(gca, 'XTickLabel', {description})
xlabel(''), ylabel('')

% Plot latents
subtightplot(2,2,3,gap,[0.1 0.25],marg_w)
h = barh([0 NaN], [latent; NaN(1,2)], 'grouped');
axis ij
xl = xlim .* [1 1.2]; xlim(xl)
dx = abs(diff(xl));
yb = [h.XOffset];
text(latent(1) + 0.01 * dx, yb(1), sprintf('\\fontsize{9}%.1f%%', explained(1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
text(latent(2) + 0.01 * dx, yb(2), sprintf('\\fontsize{9}%.1f%%', explained(2)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
title('PCA Latents'), xlabel('Variance')
set(gca, 'YTick', yb, 'YTickLabel', GLabels)
bl = gca;
bl.ActivePositionProperty = tl.ActivePositionProperty;
tlPos = get(tl, tl.ActivePositionProperty);
blPos = get(bl, bl.ActivePositionProperty);
set(bl, bl.ActivePositionProperty, [tlPos(1) blPos(2) tlPos(3) blPos(4)])

end


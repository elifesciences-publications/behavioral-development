function plotCoeffs(coeff, metricLabels, boxPCs)
%PLOTCOEFFS Plots the coefficients of a PCA analysis.
% Usage:
%   plotCoeffs(coeff, metricLabels)
%   plotCoeffs(PCs, metricLabels) % PCs is returned by pcaproj
%   plotCoeffs(..., boxPCs)
%
% See also: pcaproj, plotLatents, plotPCA, plotPC

if isstruct(coeff)
    if nargin < 2 && isfield(coeff, 'metricLabels')
        metricLabels = coeff.metricLabels;
    end
    if ~isfield(coeff, 'coeff')
        if isfield(coeff, 'pcs') && isfield(coeff.pcs, 'coeff')
            coeff = coeff.pcs;
        else
            error('Cannot find ''coeff'' subfield in input structure.')
        end
    end
    coeff = coeff.coeff;
end
if nargin < 3; boxPCs = []; end

imagesc(coeff, [-1 1]), hold on
set(gca, 'XTick', 1:size(coeff, 2))
set(gca, 'YTick', 1:size(coeff, 1))
if exist('metricLabels', 'var')
    set(gca, 'YTickLabel', metricLabels, 'TickLabelInterpreter', 'none')
end
set(gca, 'FontSize', 10)
set(gca, 'TickLength', [0,0])
xlabel('Principal Components'), ylabel('Metrics')
title('PCA Coefficients')
colorbar
colormap(redblue)

% Coefficient numbers
[x,y] = meshgrid(1:size(coeff,2), 1:size(coeff,1));
textStrings = num2str(coeff(:),'%0.2f');
textStrings = strtrim(cellstr(textStrings));
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center', 'FontSize',8);
% midValue = mean(get(gca,'CLim'));
% textColors = repmat(coeff(:) < midValue,1,3);
% set(hStrings,{'Color'},num2cell(textColors,2));

% Boxes
if ~isempty(boxPCs)
    for i = 1:numel(boxPCs)
        drawRect([boxPCs(i)-0.5 0.5 1 size(coeff,1)], 'r', 'LineWidth', 2)
    end
end

end


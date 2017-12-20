function plotPCA(pcs)
%PLOTPCA Plots the latents and the coefficients for each PC in the same figure.
% Usage:
%   plotPCA(pcs)
%
% Example:
%   plotPCA(epm.LobuleVII.Juven)
%
% See also: plotPC, plotLatents, plotCoeffs

if isfield(pcs, 'pcs')
    metricLabels = pcs.metricLabels;
    pcs = pcs.pcs;
end

subplot(1,2,1)
plotLatents(pcs)

subplot(1,2,2)
plotCoeffs(pcs, metricLabels)

% pos = get(gcf, 'Pos');
% pos(3) = pos(3) * 2;
% pos(3:4) = [1400 500];
% set(gcf, 'Pos', pos)
figsize(1400,500)

end


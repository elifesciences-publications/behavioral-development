function plotLatents(pcs)
%PLOTLATENTS Plots the latents for each PC (i.e., the variance).
% Usage:
%   plotLatents(pcs)
%
% Example:
%   plotLatents(epm.LobuleVII.Juven)
% 
% See also: pcaproj, plotCoeffs, plotPCA, plotPC

latent = [pcs.pcaLatent(:) pcs.projLatent(:)];
explained = [pcs.pcaExplained(:) pcs.projExplained(:)];
% x = 1:size(latent,1);

h = barh(latent, 'grouped');

xl = xlim .* [1 1.1]; xlim(xl)
dx = abs(diff(xl));

dy = [h.XOffset];
y = h(1).XData;
for i = 1:numel(y)
    text(latent(i,1) + 0.01 * dx, y(i) + dy(1), sprintf('%.1f%%', explained(i,1)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
    text(latent(i,2) + 0.01 * dx, y(i) + dy(2), sprintf('%.1f%%', explained(i,2)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
end

title('PCA Latents'), xlabel('Variance')
set(gca, 'YTick', y, 'YTickLabel', pcs.labels)
ylim([min(y) - 0.5, max(y) + 0.5])
axis ij
legend('Control', 'Experimental', 'Location', 'southeast')
grid on

end


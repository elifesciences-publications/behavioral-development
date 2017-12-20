function h = plotcorrgraph(Z, alpha, layout, positive_color, negative_color, rel_rho)
%PLOTCORRGRAPH Plot correlation graphs from metric tables.
% Usage:
%   plotcorrgraph(Z,...)
% 
% See also: corr_graphs

if nargin < 2 || isempty(alpha); alpha = 0.5; end
if nargin < 3 || isempty(layout); layout = 'layered'; end
if nargin < 4 || isempty(positive_color); positive_color = [1 0 0]; end
if nargin < 5 || isempty(negative_color); negative_color = [0 0 1]; end
if nargin < 6 || isempty(rel_rho); rel_rho = []; end

% Compute correlations
[rho,p] = corr(Z.Metrics,'type','Spearman','rows','pairwise');

% Offset
if ~isempty(rel_rho); rho = rho - rel_rho; end

% Create graph
G = graph(rho,'upper','OmitSelfLoops');

% Remove missing
G = G.rmedge(find(isnan(G.Edges.Weight)));

% Pull out edge data
w0 = G.Edges.Weight;
pG0 = p(sub2ind(size(p),G.Edges.EndNodes(:,1),G.Edges.EndNodes(:,2)));

% Threshold
% G = G.rmedge(find(abs(w0) < min_rho_edge));
G = G.rmedge(find(pG0 > alpha));

% Update edge data
w = G.Edges.Weight;
pG = p(sub2ind(size(p),G.Edges.EndNodes(:,1),G.Edges.EndNodes(:,2)));

% Plot
figure,figclosekey
h = G.plot('linewidth',norm2unit(abs(w)) .* 5 + 1,... % width proportional to rho
    'edgecolor',positive_color .* (w < 0) + negative_color .* (w > 0),... % color = sign
    'nodelabel',Z.Properties.VariableNames,... % node names = metrics
    'edgelabel', af(@(x)num2str(x,'%.3f'), w), ... % edge labels = rho
    'layout',layout); % node layout style
noax
axis tight
figsize(1000,800)
painters

if nargout < 1; clear h; end
end

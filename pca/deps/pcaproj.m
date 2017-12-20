function PCs = pcaproj(X, pcaIdx, projIdx, keepPCs)
%PCAPROJ Compute PCA on a subset of the data and project the rest onto the resulting PCs.
% Usage:
%   PCs = pcaproj(X, pcaIdx, projIdx)
%   PCs = pcaproj(X, pcaIdx, projIdx, keepPCs) % index of PCs to keep

if nargin < 4 || isempty(keepPCs); keepPCs = 1:size(X,2); end

% Separate subset of data used to PCA
Xpca = X(pcaIdx, :);

% Get statistics of the subset that will be used to center
mu = mean(Xpca);
sigma = std(Xpca);

% Center the data with respect to the subset of data being PCA'd
Xpca = bsxfun(@rdivide, bsxfun(@minus, Xpca, mu), sigma);
X = bsxfun(@rdivide, bsxfun(@minus, X, mu), sigma);

% PCA
[coeff, pcaScore, pcaLatent, pcaTsquared, pcaExplained] = pca(Xpca, 'Centered', false);

% Identities
assert(norm(pcaScore - Xpca * coeff) < 1e-10)
assert(norm(pcaLatent - sum(pcaScore .^ 2)' ./ size(pcaScore, 1)) < 1e-10)
assert(norm(pcaExplained - (pcaLatent ./ sum(pcaLatent)) .* 100) < 1e-10)
% assert(norm(pcaTsquared - mahal(pcaScore, pcaScore)) < 1e-10)

% Keep (some) PCs
coeff = coeff(:, keepPCs);

% Reproject all of the data onto the remaining PCs
score = X * coeff;
pcaScore = score(pcaIdx, :);   % == X(pcaIdx,:) * coeff;
projScore = score(projIdx, :); % == X(projIdx,:) * coeff;

% Compute the component variances for each subset of the data
latent = sum(score .^ 2)' ./ size(score, 1);
pcaLatent = sum(pcaScore .^ 2)' ./ size(pcaScore, 1);
projLatent = sum(projScore .^ 2)' ./ size(projScore, 1);

% Compute the fraction of explained variance by each component
explained = (latent ./ sum(latent)) * 100;
pcaExplained = (pcaLatent ./ sum(pcaLatent)) * 100;
projExplained = (projLatent ./ sum(projLatent)) * 100;

% Compute Hotelling's t-squared statistic
tsquared = mahal(score, score);
if size(pcaScore,1) > size(pcaScore,2)
    pcaTsquared = mahal(pcaScore, pcaScore);
end
if size(projScore,1) > size(projScore,2)
    projTsquared = mahal(projScore, projScore);
end

% Compute reconstruction error
Y = score * coeff';
% Y = bsxfun(@plus, bsxfun(@times, Y, sigma), mu);
[~,epsilon] = rownorm2(X - Y);

% Convenience
labels = strcat({'PC'}, strsplit(num2str(1:size(score,2))));

% Return everything in a structure
PCs = save2struct();

end


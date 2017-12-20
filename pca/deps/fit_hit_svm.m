function svm = fit_hit_svm(data, grp, polyorder)

svm = struct();

mdl = fitcsvm(data,grp,'KernelFunction','polynomial','PolynomialOrder',polyorder);

d = 0.005; % Step size of the grid
x1gv = min(mdl.X(:,1)):d:max(mdl.X(:,1));
x2gv = min(mdl.X(:,2)):d:max(mdl.X(:,2));
[x1Grid,x2Grid] = meshgrid(x1gv,x2gv);
xGrid = [x1Grid(:),x2Grid(:)];        % The grid
[~,scores] = predict(mdl,xGrid); % The scores
scores = reshape(scores(:,2),size(x1Grid));
C = contourcs(x1gv,x2gv,scores,[0 0]);

svm.C = C(1);
svm.Cpts = [svm.C.X(:) svm.C.Y(:)];
svm.mdl = mdl;
svm.scores = scores;
svm.data = data;
svm.grp = grp;

end
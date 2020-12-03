function [V, d, X_proj] = PCA(X, k)
%PCA Principal Component Analysis (PCA).
%
% [V, d, X_proj] = PCA(X, k) computes the first k eigenvectors and eigenvalues
% of the matrix X'*X/n where n is the number of rows in X.
%
% X is the dataset (N x D)
% k is the number of components to estimate or keep
%
% V is a matrix [v_1, ..., v_k] where v_i is the i-th eigenvector (D x k)
% d is the vector of the first k eigenvalues
% X_proj is the projection of X on the linear space spanned by the
% eigenvectors in V (N x k)

[n, nd] = size(X);

% empty or unexpected values for k
maxEigenVals = nd; % max number of components
if nargin<2
    k = maxEigenVals;
end
nEigs = min(k, maxEigenVals);

% main PCA
[V, D] = eigs(X'*X/n, nEigs);
d = diag(D);
d = d.*(d>0);

% redundant: eigs does eigenvalue sorting; keeping in case eig is used (for k=d);
[d, I] = sort(d, 'descend');
V = V(:,I);

% projection of data on eigenvector subspace
X_proj = X*V;

end



function [X, Y] = MixGauss(means, sigmas, n)
% Generate a mixture of p isotropic Gaussians with  diagonal covariances
%
% means: (size d x p) should be of the form [m1, ..., mp] (each mi is
% d-dimensional
% sigmas: (size p x 1) should be of the form [sigma_1; ...; sigma_p]
% n: number of points per class
%
% X: data matrix (size 2n x d)
% Y: label vector (size 2n)
%
% EXAMPLE:
% generate and visualize a 2D dataset with two classes, a 1000 points each,
% the first centered in (0,0) with variance 0.5, the second centered in
% (1,1) with variance 0.25.
%
% [X, Y] = MixGauss([[0;0],[1;1]], [0.5,0.25], 1000);
% scatter(X(:,1), X(:,2), 25, Y)
%
% see also AnisotropicMixGauss.m

[d, nComp] = size(means);

X = []; Y = [];
for c = 1:nComp
    X = [X; bsxfun(@plus, sigmas(c)*randn(d, n), means(:,c))'];
    Y = [Y; c*ones(n, 1)];
end

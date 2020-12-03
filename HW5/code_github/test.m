% [X] = generate_data('helix', 2000);
data=load('FaceData.mat');
% data=load('ScurveData.mat');
X = data.data3;
% noisestd = 0.5;
% X = X + noisestd*randn(size(X)); 
figure, scatter3(X(:,1), X(:,2), X(:,3), 5); 
title('Original dataset'); drawnow 
% no_dims = round(intrinsic_dim(X, 'MLE'));
no_dims = 3;
t=100;
sigma=0.2;
disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]); 
% [mappedX, mapping] = compute_mapping(X, 'DiffusionMaps', no_dims);
% [mappedX, mapping] = compute_mapping(X, 'PCA', no_dims);
mappedX = diffusion_maps(X, no_dims, t, sigma)
figure, scatter3(mappedX(:,1), mappedX(:,2), mappedX(:,3), 5); 
title('Result of PCA'); 
[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7); 
figure, scatter3(mappedX(:,1), mappedX(:,2),mappedX(:,3), 5); 
title('DiffusionMaps'); drawnow
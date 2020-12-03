%% PCA method for the dataset1
data_scurve=load('ScurveData.mat');
X_scurve = data_scurve.data3;
fsz = 16;
[n,d] = size(X_scurve);
figure(1);
scatter3(X_scurve(:,1),X_scurve(:,2),X_scurve(:,3),'b','filled');
title('The Scurve');

figure(2);
rng('default') % for reproducibility
[coeff,score,latent,tsquared,explained] = pca(X_scurve);
scatter(score(:,1),score(:,2),'b','+');
title('PCA');
axis equal

%% PCA method for the dataset2
data_scurve=load('ScurveData.mat');
X_scurve = data_scurve.data3;
noisestd = 0.5;
X_scurve_noised = X_scurve + noisestd*randn(size(X_scurve)); % perturb by Gaussian noise

[n,d] = size(X_scurve);
figure(3);
scatter3(X_scurve(:,1),X_scurve(:,2),X_scurve(:,3),'b','filled');
hold on
scatter3(X_scurve_noised(:,1),X_scurve_noised(:,2),X_scurve_noised(:,3),'filled');
title('The Scurve and noised Scurved with noisestd= 0.5');
legend('Scurve','noised Scurved');

figure(4);
rng('default') % for reproducibility
[coeff,score,latent,tsquared,explained] = pca(X_scurve_noised);
scatter(score(:,1),score(:,2),'b','+');
title('PCA');
axis equal

%% PCA method for the dataset3
data_facedata=load('FaceData.mat');
X_facedata = data_facedata.data3;
fsz = 16;
[n,d] = size(X_facedata);

figure(5);
rng('default') % for reproducibility
[coeff,score,latent,tsquared,explained] = pca(X_facedata);
scatter(score(:,1),score(:,2),'b','+');
title('PCA');
axis equal
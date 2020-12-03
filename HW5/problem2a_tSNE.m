%% t-SNE for unbiased scurve
data=load('ScurveData.mat');
X = data.data3;
rng('default') % for reproducibility
Y = tsne(X,'Algorithm','exact','Perplexity',10,'Distance','euclidean');
figure(1);
scatter(Y(:,1),Y(:,2),12,'b','+');
title('t-SNE');

%% t-SNE for noised scurve
data=load('ScurveData.mat');
X = data.data3;
noisestd = 0.3;
X = X + noisestd*randn(size(X)); % perturb by Gaussian noise
rng('default') % for reproducibility
Y = tsne(X,'Algorithm','exact','Perplexity',10,'Distance','euclidean');
figure(2);
scatter(Y(:,1),Y(:,2),'b','+');
title('t-SNE');

%% t-SNE for facedata
data=load('FaceData.mat');
X = data.data3;
rng('default') % for reproducibility
Y = tsne(X,'Algorithm','exact','Perplexity',10,'Distance','euclidean');
figure(3);
scatter(Y(:,1),Y(:,2),'b','+');
title('t-SNE');
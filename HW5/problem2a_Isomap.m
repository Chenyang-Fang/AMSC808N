%% Scurve without noise
% data=load('FaceData.mat');
data=load('ScurveData.mat');
X = data.data3;

fsz = 16;
% isomap method
[n,d] = size(X);
% compute pairwise distances
d = zeros(n);
e = ones(n,1);
for i = 1 : n
    d(i,:) = sqrt(sum((X - e*X(i,:)).^2,2));
end
% k-isomap
% STEP 1: find k nearest neighbors and define weighted directed graph
k = 12; % the number of nearest neighbors for computing distances
% for each point, find k nearest neighbors
ineib = zeros(n,k);
dneib = zeros(n,k);
for i = 1 : n
    [dsort,isort] = sort(d(i,:),'ascend');
    dneib(i,:) = dsort(1:k);
    ineib(i,:) = isort(1:k);
end
figure(1);
hold on;
plot3(X(:,1),X(:,2),X(:,3),'.','Markersize',15,'color','b');
daspect([1,1,1])

for i = 1 : n
    for j = 1 : k
        edge = X([i,ineib(i,j)],:);
        plot3(edge(:,1),edge(:,2),edge(:,3),'k','Linewidth',0.25);
    end
end
set(gca,'Fontsize',fsz);
view(3);


% STEP 2: compute shortest paths in the graph
D = zeros(n);
ee = ones(1,k);
g = ineib';
g = g(:)';
w = dneib';
w = w(:)';
G = sparse(kron((1:n),ee),g,w);
m = randi(n);
mf = randi(n);
c = zeros(n,3);
for i = 1 : n
    [dist,path,~] = graphshortestpath(G,i);
    D(i,:) = dist;
    if i == m
        figure(2)
        hold on
        dmax = max(dist);
        N = 1000;
        col = parula(N);
        for ii = 1 : n
            c(ii,:) = col(getcolor(dist(ii),dmax,N),:);
            plot3(X(ii,1),X(ii,2),X(ii,3),'.','Markersize',15,'color',c(ii,:));
        end
        p = path{[mf]};
        for j = 2 : length(p)
            I = [p(j-1),p(j)];
            plot3(X(I,1),X(I,2),X(I,3),'Linewidth',2,'color','r');
        end
        view(3)
        daspect([1,1,1])
        set(gca,'Fontsize',fsz);
    end
end

% STEP 3: do MDS
% symmetrize D
D = 0.5*(D + D');
Y = cmdscale(D,2);
figure(3);
hold on
for ii = 1 : n
    plot(Y(ii,1),Y(ii,2),'.','Markersize',15,'color',c(ii,:));
end
% plot edges
for i = 1 : n
    for j = 1 : k
        edge = Y([i,ineib(i,j)],:);
        plot(edge(:,1),edge(:,2),'k','Linewidth',0.25);
    end
end
% plot path
for j = 2 : length(p)
    I = [p(j-1),p(j)];
    plot(Y(I,1),Y(I,2),'Linewidth',2,'color','r');
end
set(gca,'Fontsize',fsz);
daspect([1,1,1]);

%% Scurve without noise
% data=load('FaceData.mat');
data=load('ScurveData.mat');
X = data.data3;
noisestd = 0.5;
X = X + noisestd*randn(size(X)); % perturb by Gaussian noise
fsz = 16;
% isomap method
[n,d] = size(X);
% compute pairwise distances
d = zeros(n);
e = ones(n,1);
for i = 1 : n
    d(i,:) = sqrt(sum((X - e*X(i,:)).^2,2));
end
% k-isomap
% STEP 1: find k nearest neighbors and define weighted directed graph
k = 15; % the number of nearest neighbors for computing distances
% for each point, find k nearest neighbors
ineib = zeros(n,k);
dneib = zeros(n,k);
for i = 1 : n
    [dsort,isort] = sort(d(i,:),'ascend');
    dneib(i,:) = dsort(1:k);
    ineib(i,:) = isort(1:k);
end
figure(4);
hold on;
plot3(X(:,1),X(:,2),X(:,3),'.','Markersize',15,'color','b');
daspect([1,1,1])

for i = 1 : n
    for j = 1 : k
        edge = X([i,ineib(i,j)],:);
        plot3(edge(:,1),edge(:,2),edge(:,3),'k','Linewidth',0.25);
    end
end
set(gca,'Fontsize',fsz);
view(3);


% STEP 2: compute shortest paths in the graph
D = zeros(n);
ee = ones(1,k);
g = ineib';
g = g(:)';
w = dneib';
w = w(:)';
G = sparse(kron((1:n),ee),g,w);
m = randi(n);
mf = randi(n);
c = zeros(n,3);
for i = 1 : n
    [dist,path,~] = graphshortestpath(G,i);
    D(i,:) = dist;
    if i == m
        figure(5)
        hold on
        dmax = max(dist);
        N = 1000;
        col = parula(N);
        for ii = 1 : n
            c(ii,:) = col(getcolor(dist(ii),dmax,N),:);
            plot3(X(ii,1),X(ii,2),X(ii,3),'.','Markersize',15,'color',c(ii,:));
        end
        p = path{[mf]};
        for j = 2 : length(p)
            I = [p(j-1),p(j)];
            plot3(X(I,1),X(I,2),X(I,3),'Linewidth',2,'color','r');
        end
        view(3)
        daspect([1,1,1])
        set(gca,'Fontsize',fsz);
    end
end

% STEP 3: do MDS
% symmetrize D
D = 0.5*(D + D');
Y = cmdscale(D,2);
figure(6);
hold on
for ii = 1 : n
    plot(Y(ii,1),Y(ii,2),'.','Markersize',15,'color',c(ii,:));
end
% plot edges
for i = 1 : n
    for j = 1 : k
        edge = Y([i,ineib(i,j)],:);
        plot(edge(:,1),edge(:,2),'k','Linewidth',0.25);
    end
end
% plot path
for j = 2 : length(p)
    I = [p(j-1),p(j)];
    plot(Y(I,1),Y(I,2),'Linewidth',2,'color','r');
end
set(gca,'Fontsize',fsz);
daspect([1,1,1]);
title('Isomap');
%%
function c = getcolor(u,umax,N)
c = max(1,round(N*(u/umax)));
end

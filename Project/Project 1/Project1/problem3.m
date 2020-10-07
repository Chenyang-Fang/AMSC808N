function problem3()
close all
%% read data
A2012 = readmatrix('A2012.csv');
A2016 = readmatrix('A2016.csv');
% Format for A2012 and A2016:
% FIPS, County, #DEM, #GOP, then <str> up to Unemployment Rate
str = ["Median Income", "Migration Rate", "Birth Rate",...
"Death Rate", "Bachelor Rate", "Unemployment Rate","log(#Votes)"];
%
% remove column county that is read by matlab as NaN
A2012(:,2) = [];
A2016(:,2) = [];
%% Remove rows with missing data
A = A2016;
% remove all rows with missing data
ind = find(~isfinite(A(:,2)) |  ~isfinite(A(:,3)) | ~isfinite(A(:,4)) ...
    | ~isfinite(A(:,5)) | ~isfinite(A(:,6)) | ~isfinite(A(:,7)) ...
    | ~isfinite(A(:,8)) | ~isfinite(A(:,9)));
A(ind,:) = [];
%% select CA, OR, WA, NJ, NY counties
% ind = find((A(:,1)>=6000 & A(:,1)<=6999)); % ...  %CA
%  | (A(:,1)>=53000 & A(:,1)<=53999) ...        %WA
%  | (A(:,1)>=34000 & A(:,1)<=34999) ...        %NJ  
%  | (A(:,1)>=36000 & A(:,1)<=36999) ...        %NY
%  | (A(:,1)>=41000 & A(:,1)<=41999));          %OR
% A = A(ind,:);

[n,dim] = size(A);

%% assign labels: -1 = dem, 1 = GOP
idem = find(A(:,2) >= A(:,3));
igop = find(A(:,2) < A(:,3));
num = A(:,2)+A(:,3);
label = zeros(n,1);
label(idem) = -1;
label(igop) = 1;

%% select max subset of data with equal numbers of dem and gop counties
ngop = length(igop);
ndem = length(idem);
if ngop > ndem
    rgop = randperm(ngop,ndem);
    Adem = A(idem,:);
    Agop = A(igop(rgop),:);
    A = [Adem;Agop];
else
    rdem = randperm(ndem,ngop);
    Agop = A(igop,:);
    Adem = A(idem(rdem),:);
    A = [Adem;Agop];
end  
[n,dim] = size(A)
idem = find(A(:,2) >= A(:,3));
igop = find(A(:,2) < A(:,3));
num = A(:,2)+A(:,3);
label = zeros(n,1);
label(idem) = -1;
label(igop) = 1;

%% set up data matrix and visualize
% close all
% figure;
% hold on; grid;
X = [A(:,4:9),log(num)];
X(:,1) = X(:,1)/1e4;
% select three data types that distinguish dem and gop counties the most
i1 = 1; % Median Income
i2 = 7; % log(# votes)
i3 = 5; % Bachelor Rate
% plot3(X(idem,i1),X(idem,i2),X(idem,i3),'.','color','b','Markersize',20);
% plot3(X(igop,i1),X(igop,i2),X(igop,i3),'.','color','r','Markersize',20);
% view(3)
% fsz = 16;
% set(gca,'Fontsize',fsz);
% xlabel(str(i1),'Fontsize',fsz);
% ylabel(str(i2),'Fontsize',fsz);
% zlabel(str(i3),'Fontsize',fsz);
%% rescale data to [0,1] and visualize
% figure;
% hold on; grid;
XX = X(:,[i1,i2,i3]); % data matrix
% rescale all data to [0,1]
xmin = min(XX(:,1)); xmax = max(XX(:,1));
ymin = min(XX(:,2)); ymax = max(XX(:,2));
zmin = min(XX(:,3)); zmax = max(XX(:,3));
X1 = (XX(:,1)-xmin)/(xmax-xmin);
X2 = (XX(:,2)-ymin)/(ymax-ymin);
X3 = (XX(:,3)-zmin)/(zmax-zmin);
XX = [X1,X2,X3];
% plot3(XX(idem,1),XX(idem,2),XX(idem,3),'.','color','b','Markersize',20);
% plot3(XX(igop,1),XX(igop,2),XX(igop,3),'.','color','r','Markersize',20);
% view(3)
% fsz = 16;
% set(gca,'Fontsize',fsz);
% xlabel(str(i1),'Fontsize',fsz);
% ylabel(str(i2),'Fontsize',fsz);
% zlabel(str(i3),'Fontsize',fsz);
%% set up optimization problem
[n,dim] = size(XX);
lam = 0.01;
Y = (label*ones(1,dim + 1)).*[XX,ones(n,1)];
w = [-1;-1;1;1];
fun = @(I,Y,w)fun0(I,Y,w,lam);
gfun = @(I,Y,w)gfun0(I,Y,w,lam);
Hvec = @(I,Y,w,v)Hvec0(I,Y,w,v,lam);

iter = 1000; % number of experiments for each batch size
maxiter = 2e3;
batch_set = [32 64 128];
num_batch = size(batch_set,2);
T = zeros(num_batch,2); %define the running time vector
F_SG = zeros(maxiter,num_batch);
Ave_F_SG = zeros(maxiter,num_batch);
F_SINewton = zeros(maxiter+1,num_batch);
Ave_F_SINewton = zeros(maxiter+1,num_batch);

for m = 1: num_batch
    tic
    for j = 1:iter
        w = [-1;-1;1;1];
        bsz = batch_set(1,m);
        [w,f,normf, normgrad] = SG(fun,gfun,Hvec,Y,w,bsz,maxiter);
        F_SG(:,j) = normf;
    end
    T(m,1) = toc;
    Ave_F_SG(:,m) = mean(F_SG,2);
end

figure;
hold on 
plot(Ave_F_SG(:,1));
plot(Ave_F_SG(:,2));
plot(Ave_F_SG(:,3));
xlabel('Number of iteration');
ylabel('Average function value');
legend('batch size=32','batch size=64','batch size=128');
title({['Plot of average function value f vs iteration number']
    ['for SG method ']});

%%
for m = 1: num_batch
    tic
    
    for j = 1:iter
        w = [-1;-1;1;1];
        bsz = batch_set(1,m);
        [w,f,gnorm] = SINewton2(fun,gfun,Hvec,Y,w,bsz,maxiter);
        F_SINewton(:,j) = f;
    end
    T(m,2) = toc;
    Ave_F_SINewton(:,m) = mean(F_SINewton,2);
end

figure;
hold on 
plot(Ave_F_SINewton(:,1));
plot(Ave_F_SINewton(:,2));
plot(Ave_F_SINewton(:,3));
xlabel('Number of iteration');
ylabel('Average function value');
legend('batch size=32','batch size=64','batch size=128');
title({['Plot of average function value f vs iteration number']
    ['for Subsampled inexact Newton method']});

%
figure;
hold on 
n = [32 64 128];
plot(n,T(:,1),n,T(:,2))
xlabel('batch size');
ylabel('Running time(s)')
legend('SG method','SINewton method');
title('Comparison of running time bewteen two methods');
%
figure;
hold on 
plot(Ave_F_SINewton(1:2000,1));
plot(Ave_F_SG(1:2000,1));
xlabel('Number of iteration');
ylabel('Average function value');
legend('Subsampled inexact Newton method','SG method');
title({['Comparison of average function value f vs iteration number']
    ['for two methods when batch size is 32']});

figure;
hold on 
plot(Ave_F_SINewton(1:2000,2));
plot(Ave_F_SG(1:2000,2));
xlabel('Number of iteration');
ylabel('Average function value');
legend('Subsampled inexact Newton method','SG method');
title({['Comparison of average function value f vs iteration number']
    ['for two methods when batch size is 64']});

figure;
hold on 
plot(Ave_F_SINewton(1:2000,3));
plot(Ave_F_SG(1:2000,3));
xlabel('Number of iteration');
ylabel('Average function value');
legend('Subsampled inexact Newton method','SG method');
title({['Comparison of average function value f vs iteration number']
    ['for two methods when batch size is 128']});
end
%%
function f = fun0(I,Y,w,lam)
f = sum(log(1 + exp(-Y(I,:)*w)))/length(I) + 0.5*lam*w'*w;
end
%%
function g = gfun0(I,Y,w,lam)
aux = exp(-Y(I,:)*w);
d1 = size(Y,2);
g = sum(-Y(I,:).*((aux./(1 + aux))*ones(1,d1)),1)'/length(I) + lam*w;
end
%%
function Hv = Hvec0(I,Y,w,v,lam)
aux = exp(-Y(I,:)*w);
d1 = size(Y,2);
Hv = sum(Y(I,:).*((aux.*(Y(I,:)*v)./((1+aux).^2)).*ones(1,d1)),1)' + lam*v;
end

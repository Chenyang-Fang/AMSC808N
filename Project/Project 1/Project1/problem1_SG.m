function problem1_SG()
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
ind = find((A(:,1)>=6000 & A(:,1)<=6999)); ...  %CA
%    | (A(:,1)>=53000 & A(:,1)<=53999) ...        %WA
%    | (A(:,1)>=34000 & A(:,1)<=34999) ...        %NJ  
%    | (A(:,1)>=36000 & A(:,1)<=36999) ...        %NY
%    | (A(:,1)>=41000 & A(:,1)<=41999));          %OR
A = A(ind,:);

[n,dim] = size(A);

%% assign labels: -1 = dem, 1 = GOP
idem = find(A(:,2) >= A(:,3));
igop = find(A(:,2) < A(:,3));
num = A(:,2)+A(:,3);
label = zeros(n,1);
label(idem) = -1;
label(igop) = 1;

%% set up data matrix and visualize
close all
figure;
hold on; grid;
X = [A(:,4:9),log(num)];
X(:,1) = X(:,1)/1e4;
% select three data types that distinguish dem and gop counties the most
i1 = 1; % Median Income
i2 = 7; % log(# votes)
i3 = 5; % Bachelor Rate
plot3(X(idem,i1),X(idem,i2),X(idem,i3),'.','color','b','Markersize',20);
plot3(X(igop,i1),X(igop,i2),X(igop,i3),'.','color','r','Markersize',20);
view(3)
fsz = 16;
set(gca,'Fontsize',fsz);
xlabel(str(i1),'Fontsize',fsz);
ylabel(str(i2),'Fontsize',fsz);
zlabel(str(i3),'Fontsize',fsz);

%% rescale data to [0,1] and visualize
figure;
hold on; grid;
XX = X(:,[i1,i2,i3]); % data matrix
% rescale all data to [0,1]
xmin = min(XX(:,1)); xmax = max(XX(:,1));
ymin = min(XX(:,2)); ymax = max(XX(:,2));
zmin = min(XX(:,3)); zmax = max(XX(:,3));
X1 = (XX(:,1)-xmin)/(xmax-xmin);
X2 = (XX(:,2)-ymin)/(ymax-ymin);
X3 = (XX(:,3)-zmin)/(zmax-zmin);
XX = [X1,X2,X3];
plot3(XX(idem,1),XX(idem,2),XX(idem,3),'.','color','b','Markersize',20);
plot3(XX(igop,1),XX(igop,2),XX(igop,3),'.','color','r','Markersize',20);
view(3)
fsz = 16;
set(gca,'Fontsize',fsz);
xlabel(str(i1),'Fontsize',fsz);
ylabel(str(i2),'Fontsize',fsz);
zlabel(str(i3),'Fontsize',fsz);
title({
    ['Separation of two sets of data' ] 
    ['by using SVN soft margin method'] 
    ['when batch size is 64 and n=58']
    },'Fontsize',fsz);
%% set up figure
% figure;
% hold on;
% iminus = find(label == -1);
% plot(X(iminus,1),X(iminus,2),'Linestyle','none','Marker','s','color','k');
% iplus = setdiff((1:l)',iminus);
% plot(X(iplus,1),X(iplus,2),'Linestyle','none','Marker','<','color','b');
% set(gca,'Fontsize',fsz);
% xlabel('x_1','Fontsize',fsz)
% ylabel('x_2','Fontsize',fsz)
% Lx = max(X(:,1))-min(X(:,1));
% Ly = max(X(:,2))-min(X(:,2));
% px = 0.1*Lx + 0.1;
% py = 0.1*Ly + 0.1;
% axis([min(X(:,1))-px,max(X(:,1))+px,min(X(:,2))-py,max(X(:,2))+py]);
% daspect([1,1,1]);
%% set up quadratic programming problems
[n,dim] = size(XX);
C = 1e3;
D = eye(dim+(n+1));
D(dim+1:end,dim+1:end)=0;
Y = [(label*ones(1,dim + 1)).*[XX,ones(n,dim-2)] eye(n,n);zeros(n,dim+1) eye(n,n)];
r = [ones(n,1);zeros(n,1)];
fun = @(a)(1/2)*a'*D*a + C*a;
gfun = @(a)D*a + C;
Hfun = @(a)D;

%% other format of the quadratic programming problems
% [n,dim] = size(XX);
% C = 1e3;
% e = ones(n,1);
% D = eye(dim+1); D(end,end)=0;
% DD = zeros(dim+n+1);
% DD(1:dim+1,1:dim+1)=D;
% Y = [(label*ones(1,dim + 1)).*[XX,ones(n,1)], eye(n);zeros(size(XX)+[0,1]), eye(n)];
% r = [ones(n,1);zeros(n,1)];
% fun = @(a)(1/2)*a(1:dim+1)'*D*a(1:dim+1)+C*sum(a(dim+2:end));
% gfun = @(a)[D*a(1:dim+1);C*e];
% Hfun = @(a)DD;
%% minimization using the active set method (Nocedal & Wright, section 16.5)
wb = ones(dim+n+1,1);

tic
[wb,f,~] = FindInitGuess(wb,Y,r);
if f > 1e-15
    return
end
W = find((Y*wb-r)<1e-14); % the initial active set
[wbiter, lm] = ASM(wb,gfun,Hfun,Y,r,W);
wb = wbiter(:,end);
toc

fprintf('w=[%d,%d,%d],b = %d\n',wb(1),wb(2),wb(3),wb(4));

xmin = min(XX(:,1)); xmax = max(XX(:,1));
ymin = min(XX(:,2)); ymax = max(XX(:,2));
zmin = min(XX(:,3)); zmax = max(XX(:,3));
nn = 50;
[xx,yy,zz] = meshgrid(linspace(xmin,xmax,nn),linspace(ymin,ymax,nn),...
    linspace(zmin,zmax,nn));
plane = wb(1)*xx+wb(2)*yy+wb(3)*zz+wb(4);
p = patch(isosurface(xx,yy,zz,plane,0));
p.FaceColor = 'green';
p.EdgeColor = 'none';
camlight 
lighting gouraud
alpha(0.3);

%% plot the running time
figure;
T1 = [0.206874 0.183414 0.217364 0.218842 0.217578];
T2 = [4.708526 14.035927 18.27669 58.414239 82.120475];
n = [58,97,118,180,216];
loglog(n,T1,n,T2);
xlabel('Number of counties n');
ylabel('Running time(s)');
legend('Loss function','SVN with soft margin');
title('Comparison of running time between two methods in loglog plot')

%%
% figure;
% hold on;
% grid;
% niter = length(f);
% plot((0:niter-1)',f,'Linewidth',2);
% set(gca,'Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('f','Fontsize',fsz);
%%
% figure;
% hold on;
% grid;
% niter = length(gnorm);
% plot((0:niter-1)',gnorm,'Linewidth',2);
% set(gca,'Fontsize',fsz);
% set(gca,'YScale','log');
% xlabel('k','Fontsize',fsz);
% ylabel('|| stoch grad f||','Fontsize',fsz);
end
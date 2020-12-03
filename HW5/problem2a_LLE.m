%% LLE for unbiased scurve
% GENERATE SAMPLED DATA
  data=load('ScurveData.mat');
  X = data.data3;
  [N,D] = size(X);
  K=12;
  d=2; 

  
clf; colordef none; colormap jet; set(gcf,'Position',[200,400,620,200]);

% PLOT TRUE MANIFOLD
  tt = [-1:0.1:0.5]*pi; uu = tt(end:-1:1); hh = [0:0.1:1]*5;
  xx = [cos(tt) -cos(uu)]'*ones(size(hh));
  yy = ones(size([tt uu]))'*hh;
  zz = [sin(tt) 2-sin(uu)]'*ones(size(hh));
  cc = [tt uu]' * ones(size(hh));

  figure(1);
  surf(xx,yy,zz,cc);
  view([12 -20 3]); grid off; axis off; hold on;
  lnx=-1*[1,1,1;1,-1,1]; lny=[0,0,0;5,0,0]; lnz=-1*[1,1,1;1,1,-1];
  lnh=line(lnx,lny,lnz);
  set(lnh,'Color',[1,1,1],'LineWidth',2,'LineStyle','-','Clipping','off');
  axis([-1,1,0,5,-1,3]); drawnow;

% SCATTERPLOT OF SAMPLED DATA
  figure(2);
  scatter3(X(:,1),X(:,2),X(:,3),12,'+');
  view([12 -20 3]); grid off; axis off;  hold on;
  lnh=line(lnx,lny,lnz);
  set(lnh,'Color',[1,1,1],'LineWidth',2,'LineStyle','-','Clipping','off');
  axis([-1,1,0,5,-1,3]); drawnow;

% RUN LLE ALGORITHM
  Y=lle(X',K,d);

% SCATTERPLOT OF EMBEDDING
  figure(3);
  scatter(Y(1,:),Y(2,:),12,'b','+');
  grid off; 
  title('LLE');
  %% LLE for noised scurve
  % GENERATE SAMPLED DATA
  data=load('ScurveData.mat');
  X = data.data3;
  noisestd = 0.3;
  X = X + noisestd*randn(size(X)); % perturb by Gaussian noise
  [N,D] = size(X);
  K=12;
  d=2; 

  
clf; colordef none; colormap jet; set(gcf,'Position',[200,400,620,200]);

% PLOT TRUE MANIFOLD
  tt = [-1:0.1:0.5]*pi; uu = tt(end:-1:1); hh = [0:0.1:1]*5;
  xx = [cos(tt) -cos(uu)]'*ones(size(hh));
  yy = ones(size([tt uu]))'*hh;
  zz = [sin(tt) 2-sin(uu)]'*ones(size(hh));
  cc = [tt uu]' * ones(size(hh));

  figure(4);
  surf(xx,yy,zz,cc);
  view([12 -20 3]); grid off; axis off; hold on;
  lnx=-1*[1,1,1;1,-1,1]; lny=[0,0,0;5,0,0]; lnz=-1*[1,1,1;1,1,-1];
  lnh=line(lnx,lny,lnz);
  set(lnh,'Color',[1,1,1],'LineWidth',2,'LineStyle','-','Clipping','off');
  axis([-1,1,0,5,-1,3]); drawnow;

% SCATTERPLOT OF SAMPLED DATA
  figure(5);
  scatter3(X(:,1),X(:,2),X(:,3),12,'+');
  view([12 -20 3]); grid off; axis off;  hold on;
  lnh=line(lnx,lny,lnz);
  set(lnh,'Color',[1,1,1],'LineWidth',2,'LineStyle','-','Clipping','off');
  axis([-1,1,0,5,-1,3]); drawnow;

% RUN LLE ALGORITHM
  Y=lle(X',K,d);

% SCATTERPLOT OF EMBEDDING
  figure(6);
  scatter(Y(1,:),Y(2,:),12,'b','+');
  grid off; 
  title('LLE');
  
  
  %% LLE for facedata
  % GENERATE SAMPLED DATA
  data=load('FaceData.mat');
  X = data.data3;
  [N,D] = size(X);
  K=12;
  d=2; 

  
clf; colordef none; colormap jet; set(gcf,'Position',[200,400,620,200]);

% PLOT TRUE MANIFOLD
  tt = [-1:0.1:0.5]*pi; uu = tt(end:-1:1); hh = [0:0.1:1]*5;
  xx = [cos(tt) -cos(uu)]'*ones(size(hh));
  yy = ones(size([tt uu]))'*hh;
  zz = [sin(tt) 2-sin(uu)]'*ones(size(hh));
  cc = [tt uu]' * ones(size(hh));

  subplot(1,3,1); cla;
  surf(xx,yy,zz,cc);
  view([12 -20 3]); grid off; axis off; hold on;
  lnx=-1*[1,1,1;1,-1,1]; lny=[0,0,0;5,0,0]; lnz=-1*[1,1,1;1,1,-1];
  lnh=line(lnx,lny,lnz);
  set(lnh,'Color',[1,1,1],'LineWidth',2,'LineStyle','-','Clipping','off');
  axis([-1,1,0,5,-1,3]); drawnow;

% SCATTERPLOT OF SAMPLED DATA
  subplot(1,3,2); cla;
  scatter3(X(:,1),X(:,2),X(:,3),12,'+');
  view([12 -20 3]); grid off; axis off;  hold on;
  lnh=line(lnx,lny,lnz);
  set(lnh,'Color',[1,1,1],'LineWidth',2,'LineStyle','-','Clipping','off');
  axis([-1,1,0,5,-1,3]); drawnow;

% RUN LLE ALGORITHM
Y=lle(X',K,d);

% SCATTERPLOT OF EMBEDDING
  subplot(1,3,3); cla;
  scatter(Y(1,:),Y(2,:),12,'+');
  grid off; 
  set(gca,'XTick',[]); set(gca,'YTick',[]); 
function problem3()
fsz = 16; % Fontsize
nt = 10; % trial mesh is nt-by-nt
N = 10; % the number of neurons
tol = 1e-4; % stop if ||J^\top r|| <= tol
iter_max = 120;  % max number of iterations allowed
[GDf,GDg] = GD_ode(nt,N,tol,iter_max);
% [GNf,GNg] = GaussNewton_ode(nt,N,tol,iter_max);
% [LMf,LMg] = LevenbergMarquardt_ode(nt,N,tol,iter_max);


% plot(t, err);
% xlabel('x');
% ylabel('Solution accuracy');
% title('Accuracy of the computed solution')
% figure(3);clf;
% subplot(2,1,1);
% hold on;
% % plot((1:length(GDf))',GDf,'Linewidth',2,'Marker','.','Markersize',20);
% plot((1:length(GNf))',GNf,'Linewidth',2,'Marker','.','Markersize',20);
% % plot((1:length(LMf))',LMf,'Linewidth',2,'Marker','.','Markersize',20);
% legend('Gauss-Newton');
% grid;
% set(gca,'YScale','log','Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('f','Fontsize',fsz);
% subplot(2,1,2);
% hold on;
% % plot((1:length(GDg))',GDg,'Linewidth',2,'Marker','.','Markersize',20);
% plot((1:length(GNg))',GNg,'Linewidth',2,'Marker','.','Markersize',20);
% % plot((1:length(LMg))',LMg,'Linewidth',2,'Marker','.','Markersize',20);
% legend('Gradient descend','Gauss-Newton','Levenberg-Marquardt');
% grid;
% set(gca,'YScale','log','Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('|| grad f||','Fontsize',fsz);
end
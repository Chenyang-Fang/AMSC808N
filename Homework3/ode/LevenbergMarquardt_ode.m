function [fall,norg] = LevenbergMarquardt(nt,N,tol,iter_max)
fsz = 16; % fontsize
%%
Rmax = 1;
Rmin = 1e-14;
rho_good = 0.75;
rho_bad = 0.25;
eta = 0.01;
% iter_max = 10000;
% tol = 5e-3;
%% setup training mesh
% nt = 5;
t = linspace(0,1,nt+1);
xm = meshgrid(t);
I = 2:(nt+1);
xaux = xm(1,I);
xy = [xaux(:)]';
%% initial guess for parameters
% N = 10; % the number of hidden nodes
npar = 3*N;
w = ones(npar,1);
%%
[r,J] = Res_and_Jac_ode(w,xy);
f = F(r);
g = J'*r;
nor = norm(g);
R = Rmax/5; % initial trust region radius

fprintf('Initially: f = %d, nor(g) = %d\n',f,nor); 
%% The trust region BFGS method
tic

iter = 0;
flag = 1;
I = eye(length(w));
% quadratic model: m(p) = (1/2)||r||^2 + p'*J'*r + (1/2)*p'*J'*J*p;
norg = zeros(iter_max+1,0);
fall = zeros(iter_max+1,0);
norg(1) = nor;
fall(1) = f;
while nor > tol && iter < iter_max
    % solve the constrained minimization problem using dogleg strategy
    B = J'*J + (1e-12)*I;
    pstar = -B\g;
    fprintf('iter %d: ',iter);
    if norm(pstar) <= R
        p = pstar;
        fprintf('Global min of quad model\n');
    else % solve constrained minimization problem
        lam = 1;
        isub = 0;
        while 1
            B1 = B + lam*I;
            C = chol(B1);
            p = -C\(C'\g);
            np = norm(p);
            dd = abs(np - R);
            if dd < 1e-6
                break
            end
            q = C'\p;
            nq = norm(q);
            lamnew = lam +(np/nq)^2*(np - R)/R;
            if lamnew < 0
                lam = 0.5*lam;
            else
                lam = lamnew;
            end
            isub = isub + 1;
        end
        fprintf('Contraint minimization: %d substeps\n',isub);
    end
    iter = iter + 1;  
    if flag == 0
        break;
    end
    % assess the progress
    wnew = w + p;
    [rnew, Jnew] = Res_and_Jac_ode(wnew,xy);
    mnew = 0.5*r'*r + g'*p + 0.5*p'*B*p;
    fnew = F(rnew);
    rho = (f - fnew + 1e-14)/(f - mnew + 1e-14);
    

    % adjust the trust region
    if rho < rho_bad
        R = max([0.25*R,Rmin]);
    else
        if rho > rho_good && abs(norm(p) - R) < tol
            R = min([Rmax,2*R]);
        end
    end
    % accept or reject step
    if rho > eta            
        w = wnew;
        r = rnew;
        J = Jnew;
        f = fnew;
        g = J'*r;
        nor = norm(g);        
%         fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
    end
    norg(iter+1) = nor;
    fall(iter+1) = f;
end
fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
cputime = toc;
fprintf('CPUtime = %d, iter = %d\n',cputime,iter);
%% visualize the solution
nt = 201;
t = linspace(0,2,nt);
% [xm,ym] = meshgrid(t,t);
[fun,~,~,~] = ActivationFun();
[v,W,u] = param(w);
[f0,f1,~,~,h,~,~,~,exact_sol] = setup_ode();
% A = @(x,y)(1-x).*f0(y) + x.*f1(y) + (1-y).*(g0(x)-((1-x)*f0(0)+x*f1(0))) + ...
%      y.*(g1(x)-((1-x)*f0(1)+x*f1(1)));
% B = h(xm).*h(ym);
NNfun = zeros(1,nt);
for i = 1 : nt
       NNfun(i) = v'*fun(W.*t(i) + u);
end
sol = 1 + t.*NNfun;
exact_sol = @(x)exp((-x.^2)./2)./(1+ x + x.^3)+x.^2;
esol = exact_sol(t);
err = sol - esol;
fprintf('max|err| = %d, L2 err = %d\n',max(max(abs(err))),norm(err(:)));

%
% figure(1);clf;
% hold on
% contourf(t,t,sol,linspace(min(min(sol)),max(max(sol)),20));
% plot(xy(1,:),xy(2,:),'w.','Markersize',20);
% colorbar;
% set(gca,'Fontsize',fsz);
% xlabel('x','Fontsize',fsz);
% ylabel('y','Fontsize',fsz);
% 
% %
% figure(2);clf;
% hold on
% contourf(t,t,err,linspace(min(min(err)),max(max(err)),20));
% plot(xy(1,:),xy(2,:),'w.','Markersize',20);
% colorbar;
% set(gca,'Fontsize',fsz);
% xlabel('x','Fontsize',fsz);
% ylabel('y','Fontsize',fsz);
% %
% figure(3);clf;
% subplot(2,1,1);
% fall(iter+2:end) = [];
% plot((1:iter+1)',fall,'Linewidth',2,'Marker','.','Markersize',20);
% grid;
% set(gca,'YScale','log','Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('f','Fontsize',fsz);
% subplot(2,1,2);
% norg(iter+2:end) = [];
% plot((1:iter+1)',norg,'Linewidth',2,'Marker','.','Markersize',20);
% grid;
% set(gca,'YScale','log','Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('|| grad f||','Fontsize',fsz);
plot(t, err);
xlabel('x');
ylabel('Solution accuracy');
title('Accuracy of the computed solution');
title({'Accuracy of the computed solution for LevenbergMarquardt method ','and activation function is tanh'});
end

%%
function p = cauchy_point(B,g,R)
    ng = norm(g);
    ps = -g*R/ng;
    aux = g'*B*g;
    if aux <= 0
        p = ps;
    else
        p = min(ng^3/(R*aux),1);
    end
end
%%
function f = F(r)
    f = 0.5*r'*r;
end

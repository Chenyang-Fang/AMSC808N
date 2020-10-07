function gnorm = LBFGS()
%%  the Rosenbrock function and parameters
a = 5;
func = @(x,y)(1-x).^2 + a*(y - x.^2).^2;  % Rosenbrock's function
gfun = @(x)[-2*(1-x(1))-4*a*(x(2)-x(1)^2)*x(1);2*a*(x(2)-x(1)^2)]; % gradient of f
Hfun = @(x)[2 + 12*a*x(1)^2 - 4*a*x(2), -4*a*x(1); -4*a*x(1), 2*a]; % Hessian of f
gam = 0.9; % line search step factor
jmax = ceil(log(1e-14)/log(gam)); % max # of iterations in line search
eta = 0.5; % backtracking stopping criterion factor
tol = 1e-10;
m = 5; % the number of steps to keep in memory
%% 
close all
figure;
hold on; grid;
x0 = [-1.3;1.5];  %initial guess
xstar = [1;1]; % the global minimizer
[xx,yy]=meshgrid(linspace(-2,2,1000),linspace(-1.5,2,1000));
ff = func(xx,yy);
plot(xstar(1),xstar(2),'r.','Markersize',40);
daspect([1,1,1])
col = [0.4,0.2,0];
%
s = zeros(2,m);
y = zeros(2,m);
rho = zeros(1,m);
gnorm = zeros(1,1000);
%
x = x0;
g = gfun(x);
gnorm(1) = norm(g);
plot(x(1),x(2),'.','color',col,'Markersize',20);
fx = func(x(1),x(2));
contour(xx,yy,ff,[fx,fx],'k','Linewidth',1);
% first do steepest decend step
a = linesearch(x,-g,g,func,eta,gam,jmax);
xnew = x - a*g;
gnew = gfun(xnew);
s(:,1) = xnew - x;
y(:,1) = gnew - g;
rho(1) = 1/(s(:,1)'*y(:,1));
plot([x(1),xnew(1)],[x(2),xnew(2)],'Linewidth',2,'color',col);
x = xnew;
g = gnew;
nor = norm(g);
gnorm(2) = nor;
plot(x(1),x(2),'.','color',col,'Markersize',20);
fx = func(x(1),x(2));
contour(xx,yy,ff,[fx,fx],'k','Linewidth',1);
iter = 1;
while nor > tol
    if iter < m
        I = 1 : iter;
        p = finddirection(g,s(:,I),y(:,I),rho(I));
    else
        p = finddirection(g,s,y,rho);
    end
    [a,j] = linesearch(x,p,g,func,eta,gam,jmax);
    if j == jmax
        p = -g;
        [a,j] = linesearch(x,p,g,func,eta,gam,jmax);
    end
    step = a*p;
    xnew = x + step;
    plot([x(1),xnew(1)],[x(2),xnew(2)],'Linewidth',2,'color',col);
    gnew = gfun(xnew);
    s = circshift(s,[0,1]); 
    y = circshift(y,[0,1]);
    rho = circshift(rho,[0,1]);
    s(:,1) = step;
    y(:,1) = gnew - g;
    rho(1) = 1/(step'*y(:,1));
    x = xnew;
    g = gnew;
    fx = func(x(1),x(2));
    if nor > 1e-1
        contour(xx,yy,ff,[fx,fx],'k','Linewidth',1);
    end
    plot(x(1),x(2),'.','color',col,'Markersize',20);
    nor = norm(g);
    iter = iter + 1;
    gnorm(iter+1) = nor;
end
fprintf('L-BFGS: %d iterations, norm(g) = %d\n',iter,nor);
set(gca,'Fontsize',16);
xlabel('x_1','Fontsize',16);
ylabel('x_2','Fontsize',16);
gnorm(iter+1:end) = [];
end

%%
function [a,j] = linesearch(x,p,g,func,eta,gam,jmax)
    a = 1;
    f0 = func(x(1),x(2));
    aux = eta*g'*p;
    for j = 0 : jmax
        xtry = x + a*p;
        f1 = func(xtry(1),xtry(2));
        if f1 < f0 + a*aux
            break;
        else
            a = a*gam;
        end
    end
end
function [w,f,gnorm] = StochasticLBFGS(fun,gfun,Hvec,Y)
iter = 1e3;
tol = 1e-10;

n = size(Y,1);
batch_g = 64;
batch_h = 128;
m = 5; % the number of steps to keep in memory
M = 10; %# of steps of updating inverse of Hessian % batch size


w0 = [-1;-1;1;1];  %initial guess


s = zeros(4,m);
y = zeros(4,m);
rho = zeros(1,m);
gnorm = zeros(iter,1);
% f = zeros(iter,1);
%
w = w0;
I = 1:n;
g = gfun(I,Y,w);
gnorm(1) = norm(g);

a = 0.1;

% first do steepest decend step
r = randi([1,n],batch_g,1);
g = gfun(I,Y,w);

wnew = w - a.*g;
gnew = gfun(I,Y,wnew);

s(:,1) = wnew - w;
y(:,1) = gnew - g;
rho(1) = 1/(s(:,1)'*y(:,1));

w = wnew;
g = gnew;
nor = norm(g);
gnorm(1) = nor;
f = zeros(iter,1);
f(1) = fun(I,Y,w);


k = 1;

while k < iter
    stepsize = 1./(1+k);
    if k < m*M
        upbd = ceil(k./M);
        I = 1 : upbd;
        p = finddirection(g,s(:,I),y(:,I),rho(I));
    else
        p = finddirection(g,s,y,rho);
    end

    step = stepsize.*p;
    wnew = w + step;
    
    if mod(k,M) == 0
        s = circshift(s,[0,1]); 
        y = circshift(y,[0,1]);
        rho = circshift(rho,[0,1]);
        
        hess_seed1 = randi([1 n], 1, batch_h);
        hess_seed2 = randi([1 n], 1, batch_h);
        
        g = gfun(hess_seed1,Y,wnew);
        gnew = gfun(hess_seed2,Y,wnew);
        
        s(:,1) = step;
        y(:,1) = gnew - g;
        rho(1) = 1/(s(:,1)'*y(:,1));
        
    else
        grad_seed = randi([1 n], 1, batch_g);
        gnew = gfun(grad_seed,Y,w);
        
    end
    
    nor = norm(gnew);
    if nor > .5
        gnorm(k+1) = norm(g);
        f(k+1) = fun(I,Y,w);
    elseif k > 300 && nor > .1
        gnorm(k+1) = norm(g);
        f(k+1) = fun(I,Y,w);
        
    else
        w = wnew;
        g = gnew;
        gnorm(k+1) = nor;
        f(k+1) = fun(I,Y,w);
    end
    k = k+1;
end
    
end
        
        
   
%%
function [a,j] = linesearch(Y,w,p,g,fun,eta,gam,jmax)
    a = 1;
    n=size(Y);
    I = 1:n;
    f0 = fun(I,Y,w);
    aux = eta*g'*p;
    for j = 0 : jmax
        wtry = w + a*p;
        f1 = fun(I,Y,wtry);
        if f1 < f0 + a*aux
            break;
        else
            a = a*gam;
        end
    end
end
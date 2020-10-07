function [w,f,normf, normgrad] = SG(fun,gfun,Hvec,Y,w, bsz,maxiter)
alpha0 = 1;
tol = 1e-10;


n = size(Y,1);
I=1:n;
[n,~] = size(Y);

% bsz = min(n,64); % batch size
f = zeros(maxiter,1);
normf = zeros(maxiter,1);
normgrad = zeros(maxiter,1);

for i = 1 : maxiter
    I1 = randi([1,n],1,bsz);
    b = gfun(I1,Y,w);
    normgrad(i,1) = norm(b);
    alpha = alpha0./(1+i./20);
%     w = w- 1/(size(I1,1))*alpha*b;
%     w = w- alpha0*b;
%     w = w- 1/(size(I1,1))*alpha*sum(b);
    w = w- alpha.*b;
%     w = w1/(size(I1,1));
    f = fun(I,Y,w);
    normf(i,1) =norm(f);
end
end
        
        
    
    

function [f0,f1,d1f0,d1f1,h,hp,hpp,rhs,exact_sol] = setup_ode()
f0 = @(x)1;
f1 = @(x)0;
% differential operator is d^2/dx^2 + d^2/dy^2
% differential operator applied to A(x,y), the bdry term
d1f0 = @(x)0;
d1f1 = @(x)0;
% differential operator applied to B(x,y) = x(1-x)y(1-y)NN(x,y,v,W,u)
h = @(x)x;
hp = @(x)1;
hpp = 0;
% right-hand side
rhs = @(x)x.^3+2*x+x.^2*(1+3*x.^2)/(1+ x + x.^3);
exact_sol = @(x)exp((-x.^2)./2)./(1+ x + x.^3)+x.^2;
end
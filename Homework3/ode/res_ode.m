function [r,dr] = res_ode(x,v,W,u,fun,dfun,d2fun)
% BVP for the Poisson equation is setup here
%
% residual functions r and theor derivatives w.r.t. parameters dr 
% are evaluated in this function
%
% computer diff_operator(Psi(x)) - RHS(x)
% boundary functions
a = x + (1+3.*x.^2)./(1+x+x.^3);
b = x.^3 + 2.*x + x.^2.*(1+3.*x.^2)./(1+x+x.^3);
% [~,~,d1f0,d1f1,h,hp,hpp,rhs,~] = setup_ode();
% differential operator is d^2/dx^2 + d^2/dy^2
% differential operator applied to A(x,y), the bdry term
% d1A = @(x)0;
% differential operator applied to B(x,y) = x(1-x)y(1-y)NN(x,y,v,W,u)
% h1 = h(x(1));
% hp1 = hp(x(1));
[f,fx,df,dfx] = NN_ode(x,v,W,u,fun,dfun,d2fun);
% d1B = hp(h1)*f;
r = f + x*fx + a*(1+x*f) - b;
% residual r = d2A + d2B - RHS
% r = d1A(x) + d1B - rhs(x(1));
% derivative of r w.r.t. parameters
dr = (a*x+1).*df + x.*dfx;
% dr = hp*(h1)*df ;
end



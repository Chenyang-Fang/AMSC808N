function [fun,dfun,d2fun,d3fun] = ActivationFun()
% activation function and its derivatives
%% tanh
fun = @(x)tanh(x);
dfun = @(x)1./cosh(x).^2;
d2fun = @(x)-2*sinh(x)./cosh(x).^3;
d3fun = @(x)(4*sinh(x).^2-2)./cosh(x).^4;
%% sigmoid
% fun = @(x)1./(1+exp(-x));
% dfun = @(x)exp(-x)./(1+exp(-x)).^2;
% d2fun = @(x)-exp(-x)./(1+exp(-x)).^2 + 2*exp(-2*x)./(1+exp(-x)).^3;
% d3fun = @(x)exp(-x)./(1+exp(-x)).^2 - 6*exp(-2*x)./(1+exp(-x)).^3 + 6*exp(-3*x)./(1+exp(-x)).^4;
end

function [f,fx,df,dfx] = NN_ode(x,v,W,u,fun,dfun,d2fun)
%% derivatives of the network 
z = W.*x + u;
s0 = fun(z); % sigma(z)
f = v'*s0;
% W2 = W.*W;
s1 = dfun(z); % sigma'(z)
s2 = d2fun(z); % sigma''(z)
% s3 = d3fun(z); % % sigma'(z)
fx = v'*(W.*s1); % Psi_x
% fy = v'*(W(:,2).*s1); % Psi_y
% fxx = v'*(W2(:,1).*s2); % Psi_{xx}
% fyy = v'*(W2(:,2).*s2); % Psi_{yy}
%% derivatives with respect to parameters
nv1 = length(v); % nv2 must be 1
nw1 = length(W);
nu1 = length(u); % nu2 must be 1
dim = nv1 + nw1 + nu1;
df = zeros(dim,1);
dfx = zeros(dim,1);
% dfy = zeros(dim,1);
% dfxx = zeros(dim,1);
% dfyy = zeros(dim,1);
% df
df(1:nv1) = s0;
df(nv1+1 : nv1+nw1) = x.*(v.*s1); %was reshape((v.*s1)*(x'),[nw1*nw2,1])
df(nv1+nw1+1 : end) = v.*s1;
% dfx
dfx(1:nv1) = W.*s1;
dfx(nv1+1 : nv1+nw1) = x.*(v.*W.*s2)+(v.*s1);
    %reshape((v.*W(:,1).*s2)*(x') + (v.*s1)*[1,0],[nw1*nw2,1]);
dfx(nv1+nw1+1 : end) = v.*W.*s2;
% dfy
% dfy(1:nv1) = W(:,2).*s1;
% dfy(nv1+1 : nv1+nw1*nw2) = ...
%     reshape((v.*W(:,2).*s2)*(x') + (v.*s1)*[0,1],[nw1*nw2,1]);
% dfy(nv1+nw1*nw2+1 : end) = v.*W(:,2).*s2;
% % dfxx
% dfxx(1:nv1) = W2(:,1).*s2;
% dfxx(nv1+1 : nv1+nw1*nw2) = ...
%     reshape((v.*W2(:,1).*s3)*(x') + 2*(v.*W(:,1).*s2)*[1,0],[nw1*nw2,1]);
% dfxx(nv1+nw1*nw2+1 : end) = v.*W2(:,1).*s3;
% % dfyy
% dfyy(1:nv1) = W2(:,2).*s2;
% dfyy(nv1+1 : nv1+nw1*nw2) = ...
%     reshape((v.*W2(:,2).*s3)*(x') + 2*(v.*W(:,2).*s2)*[0,1],[nw1*nw2,1]);
% dfyy(nv1+nw1*nw2+1 : end) = v.*W2(:,2).*s3;
end
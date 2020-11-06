A = readmatrix('Movierankings36.csv');
[m,n] = size(A);
index = isfinite(A);
A_new = zeros(m,n);
A_new(index > 0) = A(index > 0);
k = 20;
X = rand(m,k);
Y = rand(n,k);

iter = 1000;
lambda = [2,20,50]; %test for different lambda
error_A = zeros(3,iter+1);
error_A(:,1) = norm((A_new - X*Y').*index,'fro');


for i = 1:3
    for l = 1:iter
        Xnew = updateX(A_new,index,Y,lambda(i));
        Ynew = updateY(A_new,index,X,lambda(i));
        X = Xnew;
        Y = Ynew;
        error_A(i,l+1) = norm((A_new - X*Y').*index,'fro');
    end
end
figure(1);
hold on
plot(error_A(1,:));
plot(error_A(2,:));
plot(error_A(3,:));
xlabel('number of iteration');
ylabel('Norm of objective function');
legend('lambda=0.01','lambda=0.1','lambda=1');
title('Alternative iteration method verse iteration for different lambda');

%%
M = rand(m,n);
iter = 5;
lambda = [0.01,0.1,1]; %test for different lambda
error_N = zeros(3,iter+1);
error_N(:,1) = norm((A_new - M).*index,'fro');
for i = 1:3
    for l = 1:iter
        D = M + (A_new-M).*index;
        [U,S,V] = svd(D);
        Slamb = relu(dlarray(S - lambda(i).*eye(m,n)));
        M = extractdata(U * Slamb * V');
        error_N(i,l+1) = norm((A_new - M).*index,'fro');
    end
end
figure(2);
hold on
plot(error_N(1,:));
plot(error_N(2,:));
plot(error_N(3,:));
xlabel('number of iteration');
ylabel('Norm of objective function');
legend('lambda=0.01','lambda=0.1','lambda=1');
title('Nuclear norm trick verse iteration for different lambda');

%%
function Xnew = updateX(A,Omega,Y,lambda)
[m,~] = size(A);
[~,k] = size(Y);
Xnew = zeros(m,k);
for i = 1:m
    Yi = 0.*Y;
    ai = A(i,:)';
    Yi(Omega(i,:)>0,:) = Y(Omega(i,:)>0,:);
    H = Yi' * Yi + lambda.* eye(k);
    r = Yi' * ai;
    Xnew(i,:) = H\r;
end
end

function Ynew = updateY(A,Omega,X,lambda)
[~,n] = size(A);
[~,k] = size(X);
Ynew = zeros(n,k);
for i = 1:n
    Xi = 0.*X;
    ai = A(:,i);
    Xi(Omega(:,i)>0,:) = X(Omega(:,i)>0,:);
    H = Xi' * Xi + lambda.* eye(k);
    r = Xi' * ai;
    Ynew(i,:) = H\r;
end
end
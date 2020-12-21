%% problem part1: plot the flat region
j = 0:5;
x = pi*j/10;
a = linspace(-10,10,100);
b = zeros(6,100);
min_global_x = 0.8613;
min_global_y = 0.3735;

for j = 1:5
    b(j,:) = a*j*pi/10;
end
a_ray = linspace(10/pi,10,100);
b_ray = a_ray*pi/2-1;

hold on 
plot(a,b(1,:),'--',a,b(2,:),'--',a,b(3,:),'--',a,b(4,:),'--',a,b(5,:),'--',a,b(6,:),'--');
scatter(min_global_x, min_global_y, 12,'filled','r');
plot(a_ray,b_ray,'r','LineWidth',2);
text(-2, 12, '\Omega_1');
text(4, 6.8, '\rightarrow \Omega_2');
legend('j=1','j=2','j=3','j=4','j=5','j=0','global minimizer','b = a*\pi/2-1')
xlabel('a');
ylabel('b');
title('The set of stationary point of f');

%% problem1 part 2 plot the direction on the figure
j = 0:5;
x = pi*j/10;
a = linspace(0,5,100);
b = zeros(6,100);
min_global_x = 0.8613;
min_global_y = 0.3735;

a_ini = 1;
b_ini = 0;

for j = 1:5
    b(j,:) = a*j*pi/10;
end
a_ray = linspace(10/pi,5,100);
b_ray = a_ray*pi/2-1;

g1 = [1,0.5891, 0.1782, 0];
g2 = [0,0.3949, 0.7898, 0.9611];


hold on 
plot(a,b(1,:),'b--',a,b(2,:),'b--',a,b(3,:),'b--',a,b(4,:),'b--',a,b(5,:),'b',a,b(6,:),'b--');
scatter(min_global_x, min_global_y, 12,'filled','r');
scatter(a_ini, b_ini, 12,'filled','b');
plot(a_ray,b_ray,'r','LineWidth',2);
plot(g1,g2,'black');
xlabel('a');
ylabel('b');
title('The initial steepest direction of function f(a,b) at (1,0)');

stepsize_opt = (pi/2)/(0.3949+pi/2*0.4109);

%% problem1 part 2_1 2_2 \alpha^* and 0.99*\alpha^*
j = 0:5;
x = pi*j/10;
a = 1;
b = 0;

dfda = guna(a,b,x);
dfdb = gunb(a,b,x);
J = [dfda dfdb];
S = -(J);
maxiter = 500;
fval = zeros(1,maxiter+1);
fval(1,1) = fun(a,b,x);
tol = 1e-3;
for i = 1: maxiter 
    error = abs(fval(1,i)-0.1405);
    if  error > tol
        alpha = 0.99* 1.50996; % constant stepsize 
        a = a + alpha*S(1);
        b = b + alpha*S(2);
        dfda = guna(a,b,x);
        dfdb = gunb(a,b,x);
        J = [dfda dfdb]; % Updated Gradient
        S = -(J); % New Search Direction
        fval(1,i+1)= fun(a,b,x);
    else
        fval(1,i+1) = fval(1,i);
    end
end
plot(fval);
xlabel('number of iterations');
ylabel('Value of loss function f(a,b)');
title('Loss function value verse number of iterations with 0.99*\alpha^*');

%% problem1 part 2_3 global minimizer
j = 0:5;
x = pi*j/10;
a = 1;
b = 0;

dfda = guna(a,b,x);
dfdb = gunb(a,b,x);
J = [dfda dfdb];
S = -(J);
maxiter = 500;
fval = zeros(1,maxiter+1);
fval(1,1) = fun(a,b,x);
tol = 1e-5;
a_set = zeros(1,maxiter+1);
b_set = zeros(1,maxiter+1);
a_set(1,1) = a;
b_set(1,1) = b;
for i = 1: maxiter 
    error = abs(fval(1,i)-0.1405);
    if  error > tol
        alpha = 1.3; % constant stepsize 
        a = a + alpha*S(1);
        b = b + alpha*S(2);
        a_set(1,i+1) = a;
        b_set(1,i+1) = b;
        dfda = guna(a,b,x);
        dfdb = gunb(a,b,x);
        J = [dfda dfdb]; % Updated Gradient
        S = -(J); % New Search Direction
        fval(1,i+1)= fun(a,b,x);
    else
        fval(1,i+1) = fval(1,i);
    end
end
a_base = 1:maxiter;
b_base = ones(1,maxiter)*3.6327e-4;
figure(1)
hold on 
plot(fval);
plot(a_base,b_base);
xlabel('number of iterations');
ylabel('Value of loss function f(a,b)');
legend('numerical value of loss function', 'theoretical value of loss function')
title('Loss function value verse number of iterations with \alpha=1.3');

min_global_x = 0.8613;
min_global_y = 0.3735;
figure(2)
hold on 
plot(min_global_x,min_global_y,'r+','MarkerSize',16);
plot(a_set,b_set,'b');
xlabel('a');
ylabel('b');
title('Trajectory of iterates');
%% problem 1 part 3: sgd
alpha0=1.5;
batch_size=1;
j = 0:5;
x = pi*j/10;
n = length(x);
maxiter = 3000;
normgrad = zeros(1,maxiter);
f = zeros(1,maxiter);
a=1;
b=0;
a_set = zeros(1,maxiter+1);
b_set = zeros(1,maxiter+1);
a_set(1,1) = a;
b_set(1,1) = b;
for i = 1 : maxiter
    n_randi = randi([1,n],1,batch_size);
    x_randi = x(n_randi);
    dfda = guna(a,b,x_randi);
    dfdb = gunb(a,b,x_randi);
    J = [dfda dfdb];
    normgrad(1,i) = norm(J);
    alpha = alpha0./(1+i./100);
    a = a - alpha.*J(1);
    b = b - alpha.*J(2);
    a_set(1,i+1) = a;
    b_set(1,i+1) = b;
    f(1,i) = fun(a,b,x);
end
a_base = 1:maxiter;
b_base = ones(1,maxiter)*3.6327e-4;
figure(1)
hold on 
plot(f);
plot(a_base,b_base);
xlabel('number of iterations');
ylabel('Value of loss function f(a,b)');
legend('numerical value of loss function', 'theoretical value of loss function')
title('Loss function value verse number of iterations with stepsize strategy');

min_global_x = 0.8613;
min_global_y = 0.3735;
figure(2)
hold on 
plot(min_global_x,min_global_y,'r+','MarkerSize',16);
plot(a_set,b_set,'b');
xlabel('a');
ylabel('b');
title('Trajectory of iterates');
%% defined function
function f = fun(a,b,x)
n = length(x);
g = @(x)1-cos(x);
relu = @(a,b,x)max(a*x-b,0);
t = relu(a,b,x);
f=0;
for i = 1:n
f = f + 1/12 * (t(i)-g(x(i))).^2;
end
end

function dfda = guna(a,b,x)
g = @(x)1-cos(x);
relu = @(a,b,x)max(a*x-b,0);
t = relu(a,b,x);
dfda = 0;
n = length(x);
for i = 1:n
    if t(i)>0
        dfda = dfda + 1/6 *(t(i)-g(x(i)))*x(i);
    else
        dfda = dfda + 0;
    end
end
end

function dfdb = gunb(a,b,x)
g = @(x)1-cos(x);
relu = @(a,b,x)max(a*x-b,0);
t = relu(a,b,x);
dfdb = 0;
n = length(x);
for i = 1:n
    if t(i)>0
        dfdb = dfdb + 1/6 *(t(i)-g(x(i)))*(-1);
    else
        dfdb = dfdb + 0;
    end
end
end



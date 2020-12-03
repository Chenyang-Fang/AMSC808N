%% experimental solution
n = 1000;
r = 100;
m = 20;
node = 1:n;
% z = linspace(0,4,m);
z = 4/m:4/m:4;
Smean = zeros(1,m);
smean = zeros(1,m);
Smax = zeros(1,r);
savg = zeros(1,r);
for i =1:m
   for j = 1:r
       p = z(i)/(n-1); 
       rand('seed',100); % reseed so you get a similar picture
       G = rand(n,n) < p;
       G = triu(G,1);
       G = G + G';
       A = graph(G);
%        B = dfs(A);

       [bin, binsize] = conncomp(A);
       [Smax(1,j), int_max] = max(binsize);
       idx = node(~(bin == int_max));
       v = randsample(idx, 1);
       savg(1,j) = bin(v);
%        savg(1,j) = mean(bin);
   end
   Smean(1,i) = mean(Smax(1,:));
   smean(1,i) = mean(savg(1,:));
end


%% the numerical solution
Sext = zeros(1,m);
z = 4/m:4/m:4;
for k = 1:m
    param = z(k);
    syms x
    eqnLeft = x;
    eqnRight = 1 - exp(-param*x);
    Sext(k) = vpasolve(eqnLeft == eqnRight, x, 1);
end

sext = 1./(1 - z + z.*Sext);

%% generate the figures
figure(1);
hold on
plot(z,Smean/1000);
plot(z,Sext);
title('Fraction for different size random graph G(n,p)')
xlabel('z');
ylabel('S(z)');
legend('Experimental result', 'Theoretical result');

figure(2)
loglog(z,(smean/10),z,(sext));
title('Acerage size of the non-giant component in loglog form')
xlabel('z');
ylabel('s(z)');
legend('Experimental result', 'Theoretical result');

figure(3);
hold on
plot(smean);
plot(sext);
title('Acerage size of the non-giant component')
xlabel('z');
ylabel('s(z)');
legend('Experimental result', 'Theoretical result');

%% dfs function

function E = dfs(G)
N = numnodes(G);
color = cell(N,1);
color(:) = {'white'};
color = string(color);
T = table(color);
T.parent = ones(N,1)*Inf;
T.d = zeros(N,1);
T.f = zeros(N,1);
time = 0;
for u = 1:N
    if T.color(u) == "white"
        time = dfs_visit(G,u,time);
    end
end

ind = ~isinf(T.parent);
nodes = (1:N)';
E = graph(T.parent(ind), nodes(ind));
end

function time = dfs_visit(G,u,time)
global T
time = time + 1;
T.d(u) = time;
T.color(u) = "gray";

adj_v = neighbors(G,u);
g = length(adj_v);
for k = 1:g
    v = adj_v(k);
    if T.color(v) == "white"
        T.parent(v) = u;
        time = dfs_visit(G,v,time);
    end
end

T.color(u) = "black";
time = time + 1;
T.f(u) = time;
end
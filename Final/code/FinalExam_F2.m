%% Problem 2(a) calculate the S
a=2.2;
G1 = @(x)polylog(a-1,x)./(x*polylog(a-1,1));
G0 = @(x)polylog(a,x)/polylog(a,1);
u = fzero(@(x)G1(x)-x,0.2);
S = 1-G0(u);
%% Problem 2(b) calculate the S
a=2.2;
T = 0.4;
G1 = @(x)polylog(a-1,x)./(x*polylog(a-1,1));
G0 = @(x)polylog(a,x)/polylog(a,1);
u = fzero(@(x)G1(1-T+T*x)-x,0.3);
S = 1-G0(1-T+T*u);

%% problem 3(a)
a=2.2;
n = 1e4;
maxiter = 100;
binsize_max = zeros(1,maxiter);
for i = 1:maxiter
    [G,edges,K,p] = MakePowerLawRandomGraph(n,a);
    A = graph(G);
    [bin, binsize] = conncomp(A);
    binsize_max(1,i) = max(binsize);
end
binsize_avg = mean(binsize_max)/n;

%% problem 3(b)
a=2.2;
n = 1e4;
maxiter = 100;
binsize_max = zeros(1,maxiter);
T=0.4;


for i = 1:maxiter
    a=2.2;
    n = 1e4;
    [G,edges,K,p] = MakePowerLawRandomGraph(n,a);
    n_edge = size(edges,1);
    for j = 1:n_edge
        M=binornd(1,T);
        if M == 1
            G_temp1 = edges(j,1);
            G_temp2 = edges(j,2);
            G(G_temp1,G_temp2) = 1;
            G(G_temp2,G_temp1) = 1;
        else 
            G_temp1 = edges(j,1);
            G_temp2 = edges(j,2);
            G(G_temp1,G_temp2) = 0;
            G(G_temp2,G_temp1) = 0;
        end
    end
    A = graph(G);
    [bin, binsize] = conncomp(A);
    binsize_max(1,i) = max(binsize);
end
binsize_avg = mean(binsize_max)/n;


%% %% problem 3(c)
a=2.2;
n = 1e4;
maxiter = 100;
T=0.01:0.01:0.4;
T_size = size(T,2);

binsize_max = zeros(1,maxiter);
binsize_avg = zeros(1,T_size);
for m = 1: T_size
    T_iter = T(m);
    for i = 1:maxiter
        a=2.2;
        n = 1e4;
        [G,edges,K,p] = MakePowerLawRandomGraph(n,a);
        n_edge = size(edges,1);
        for j = 1:n_edge
            M=binornd(1,T_iter);
            if M == 1
                G_temp1 = edges(j,1);
                G_temp2 = edges(j,2);
                G(G_temp1,G_temp2) = 1;
                G(G_temp2,G_temp1) = 1;
            else 
                G_temp1 = edges(j,1);
                G_temp2 = edges(j,2);
                G(G_temp1,G_temp2) = 0;
                G(G_temp2,G_temp1) = 0;
            end
        end
        A = graph(G);
        [bin, binsize] = conncomp(A);
        binsize_max(1,i) = max(binsize);
    end
    binsize_avg(1,m) = mean(binsize_max)/n;
end
x_ini = [0.01 0.03 0.4];
y_ini = [0.001222 0.007137 0.0952];
hold on 
plot(T,binsize_avg,'LineWidth',3);
plot(x_ini, y_ini);
xlabel('T');
ylabel('Average fraction of nodes affected');
legend('Actual transimission', 'Transimission speed T<0.03');
title('Average fraction of nodes affected verse transmissibility');

%% problem 4
T = 0.4;
alpha = 2.2;
n = 1e3;
[G,edges,K,p] = MakePowerLawRandomGraph(n,a);
A = graph(G);

maxiter = 100;
step = 1000;
inf = zeros(maxiter,step);

for i = 1: maxiter
node_init = randi(n);
inf_iter = zeros(1,step);
color = cell(n,1);
color(:) = {'white'};
color = string(color);
vertex = 1:n;

for k = 1:step
    for node_infected = node_init 
        color(node_infected) = 'gray';
        nbhd = neighbors(A,node_infected);
        for node = nbhd
            coin = binornd(1,T);
            if coin == 1
                color(node) = 'black'; 
            else
                continue;
            end
        end
    end
    node_black = (color == "black");
    node_gray = (color == "gray");
    inf_iter(k) = sum(node_black) + sum(node_gray);
    color(node_gray) = "white";
    node_init = vertex(node_black);
end

inf(i, :) = inf_iter;
end
mean_infection = mean(inf,1);
plot(mean_infection./n);
xlabel('Time step');
ylabel(' Average fraction of infected nodes');
title('Evolution of average fraction of affected nodes verse time');






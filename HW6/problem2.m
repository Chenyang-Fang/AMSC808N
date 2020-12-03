%% experimental solution
p = [10,11,12,13];
m = size(p);
n = 2.^p;
r = 100;
z = 4;
bin_min = zeros(1,r);
bin_min_sample = zeros(1,r);
bin_avg = zeros(1,4);
D = zeros(1,100);
for i = 1:4
    for j = 1:r
        p_val = p(i);
        n_val = n(i);
        q = z/(n_val-1);
        G = create_ER_Graph(n_val, q);
        A = graph(G);
        v = randi([1 n_val],1,r); 
%         B = bfs(A,k);
        for k = 1:r
            [path, path_dist] = shortestpath(A,1,v(k));
            if isfinite(path_dist)
                bin_min(k) = path_dist;
            end
        end
        bin_min_sample(1,j) = mean(nonzeros(bin_min));
    end
    bin_avg(1,i) = mean(bin_min_sample);
end


%% plot
lest = (10:13)*log(2)/log(4);
hold on 
plot(p,bin_avg);
plot(p, lest);
title('Shortest path length for different size random graph G(n,p)');
xlabel('n');
ylabel('average path length l(n)');
legend('Experimental result', 'Theoretical result');

%% bfs 
function [A, component] = bfs(G,s)
N = numnodes(G);
color = cell(N,1);
parent = cell(N,1);
color(:) = {'white'};
parent(:) = {'nil'};
color = string(color);
parent = string(parent);
T = table(color, parent);
T.d = ones(N,1)*Inf;
T.color(s) = "gray";
T.d(s) = 0;
Q = [];
Q(1) = s;

while ~isempty(Q)
    u = Q(1);
    adj_v = neighbors(G,u);
    g = length(adj_v);
    for k = 1:g
        v = adj_v(k);
        if T.color(v) == "white"
            T.color(v) = "gray";
            T.d(v) = T.d(u) + 1;
            T.parent(v) = u;
            Q(end+1) = v;
        end
    end
    T.color(u) = "black";
    Q(1) = [];
end
end

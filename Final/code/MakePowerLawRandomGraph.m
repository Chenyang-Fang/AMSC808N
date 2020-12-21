function [G,edges,K,p] = MakePowerLawRandomGraph(n,a)
% Input:
% n = the number of vertices
% a = the exponent for the power-law degree distribution: p_k ~ k^{-a}
% Output:
% G = the adjacency matrix
% edges = the list of edges
% K = the maximal degree
% p = the probability distribution
%
K = round(n^(1/(a-1))); % the maximal degree;
k = (1:K)';
p = k.^(-a);
psum = sum(p);
p = p/psum; % the desired probabilities
% assign degrees to vertices according to p
aux = round(cumsum(n*p));
w = zeros(n,1);
prev = 1;
for j = 1 : K
    if aux(j) >= prev
        w(prev:aux(j)) = j;
        prev = aux(j)+1;
    end
end
% w = vector of degrees
wcumsum = cumsum(w);
m = wcumsum(end); % (the total number of edges)*2
% m must be even. Make it even if it is odd
if mod(m,2) == 1
    w(end) = w(end) + 1;
    wcumsum = cumsum(w);
end
m = wcumsum(end);
% enumerate stubs and assign the to vertices according to their degrees
stubs = zeros(m,1);
prev = 1;
for j = 1 : n
    stubs(prev:wcumsum(j)) = j;
    prev = wcumsum(j) + 1;
end
% randomly match stubs
s0 = randperm(m);
s1 = s0(1:m/2);
s2 = s0(m/2+1:end);
edges = zeros(m/2,2);
edges(:,1) = stubs(s1);
edges(:,2) = stubs(s2);
% remove selfloops
ind = find(edges(:,1)==edges(:,2));
edges(ind,:) = [];
% remove repeated edges
edges = sort(edges,2);
edges = unique(edges,'rows');
% form the adjacency matrix
G = sparse(edges(:,1),edges(:,2),ones(size(edges,1),1),n,n);
G = spones(G + G');
% calculate the resulting degree distribution in G and plot its bar graph
deg = sum(G,2);
ndeg = zeros(K,1);
for j = 1 : K
    ndeg(j) = length(find(deg == j));
end
sdeg = sum(ndeg);
%
% figure(2);clf;
% hold on
% bar(k,ndeg/sdeg);
% plot(k,p,'.','Markersize',20);
% set(gca,'Yscale','log','Xscale','log','fontsize',16);
end

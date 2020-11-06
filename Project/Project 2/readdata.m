function [M,y] = readdata()
%% read data
fid=fopen('vectors.txt','r');
A = fscanf(fid,'%i\n');
ind = find(A>100000); % find entries of A with document IDs
n = length(ind); % the number of documents
la = length(A);
I = 1:la;
II = setdiff(I,ind);
d = max(A(II)); % the number of words in the dictionary
M = zeros(n,d);
y = zeros(n,1); % y = -1 => category 1; y = 1 => category 2 
% define M and y
for j = 1 : n
    i = ind(j); 
    y(j) = A(i+1); 
    if j<n 
        iend = ind(j+1)-1; 
    else
        iend = length(A);
    end
    M(j,A(i+2:2:iend-1)) = 1; A(i+3:2:iend);
end
i1 = find(y==-1);
i2 = find(y==1);
ii = find(M>0);
n1 = length(i1);
n2 = length(i2);

fprintf('Class 1: %d items\nClass 2: %d items\n',n1,n2);
fprintf('M is %d-by-%d, fraction of nonzero entries: %d\n',n,d,length(ii)/(n*d));
end

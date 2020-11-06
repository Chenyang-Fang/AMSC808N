function C = columnselect(A,V,k,c)
[m,n] = size(A);
C = zeros(m,n);
for j = 1:n
    score = sum(V(j,1:k).^2)./k; 
    p = min(1,c*score);
    sample = binornd(1,p);
    if sample == 1
        C(:,j) = A(:,j);
    end
end
end

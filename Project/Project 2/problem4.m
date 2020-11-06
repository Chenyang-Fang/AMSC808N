%%
[M,y] = readdata();
florida = find(y == -1);
indiana = find(y == 1);
k = 5;
[m,n] = size(M);
%produce a ecomonic size decomposition of matrix M
[U,S,V] = svd(M,'econ');

% normalized statistical leverage score
score = zeros(1,n);
for j = 1:n
    score(j) = sum(V(j,1:k).^2)./k; 
end
[~,index] = maxk(score,k);
fprintf('Indices with 5 maximal leverage scores:%d\n',index);

%%
ind = readmatrix('pos_selection.txt');
M_sub = M(:,ind);
[m_sub,n_sub] = size(M_sub);
[U_sub,S,V_sub] = svd(M_sub,'econ');
k = 5;

score_sub = zeros(1,n_sub);
for j = 1:n_sub
    score_sub(j) = sum(V_sub(j,1:k).^2)./k;
end

[~,index_select] = maxk(score_sub,k);
index_select_original = ind(index_select);
fprintf('Indices with 5 maximal leverage scores %d:\n',index_select_original);
%%
PCA1 = V_sub(:,1:2);
projection1 = M_sub * PCA1;

Msmall = M_sub(:,index_select);
[Usmall,Ssmall,Vsmall] = svd(Msmall,'econ');
PCA2 = Vsmall(:,1:2);
projection2 = Msmall * PCA2;

figure(1);
hold on
scatter(projection1(florida,1),projection1(florida,2));
scatter(projection1(indiana,1),projection1(indiana,2));
xlabel('v1');
ylabel('v2');
legend('Florida','Indiana');
title('Principal component analysis classification');
figure(2);
hold on
scatter(projection2(florida,1),projection2(florida,2));
scatter(projection2(indiana,1),projection2(indiana,2));
xlabel('v1');
ylabel('v2');
legend('Florida','Indiana');
title('Principal component analysis classification');
A = importdata('MovieRankings36.csv');
rows = [1,2,3,4,8,13,16,25,29];

columns = find(isfinite(A(1,:)) & isfinite(A(2,:)) & isfinite(A(3,:)) ...
    & isfinite(A(4,:)) & isfinite(A(8,:)) & isfinite(A(13,:)) ...
    & isfinite(A(16,:)) & isfinite(A(25,:)) & isfinite(A(29,:)));
rank = 15;
A = A(rows,columns);
[m,n] = size(A);
A = dlarray(A);

[Fnorm_pgd] = nmf_pgd(A, rank);
[Fnorm_ls] = nmf_ls(A, rank);
[Fnorm_hybrid] = nmf_hybrid(A, rank);
figure(1);
hold on
plot(Fnorm_pgd);
title('Frobenius norm square verse iteration number when cluster k=10');
xlabel('number of iteration');
ylabel('Frobenius norm');

figure(2);
hold on 
plot(Fnorm_ls);
title('Frobenius norm square verse iteration number when cluster k=10');
xlabel('number of iteration');
ylabel('Frobenius norm');

figure(3);
hold on 
plot(Fnorm_hybrid);
title('Frobenius norm square verse iteration number when cluster k=10');
xlabel('number of iteration');
ylabel('Frobenius norm');
% read data
[M,y] = readdata();
[U1,S1,V1] = svd(M,'econ');
[U2,S2,V2] = svd(M','econ');
sigma = diag(S1);
a = 8; k = 10;
pho1 = zeros(1,k-1);
pho2 = zeros(k-1,a);
num_run = 100;


for i = 2:k
    pho1(i-1) = sqrt(norm(sigma(i+1:end),'fro'));
    for j = 1:a
        c = j*i;
        r = j*i;
        iter = 1;
        sample = zeros(1,num_run);
        parfor iter = 1:num_run
            tic
            C = columnselect(M,V1,i,c);
            R = columnselect(M',V2,i,r)';
            U = pinv(C) * M * pinv(R);
            sample(iter) = sqrt(sum((M - C*U*R).^2,'all'));
        end
        pho2(i-1,j) = mean(sample);
    end
end

function [Fnorm] = nmf_hybrid(A, rank)
% Projected gradient descend for non-negative matrix factorization
%
% The problem of interest is defined as 
%    min || A- WH ||_F^2,
%    where {A,W,H} >0.
%
% Given a non-negative matrix A, factorized non-negative matrice
%
% Inputs:
%    A    : (m x n) non-negative matrix to factorize
%    rank : rank
%
% Outputs:

% set dimensions and samples
[m,n] = size(A);
% initialize factors 
W = rand(m,rank);
H = rand(rank,n);
% initialize
Fnorm = zeros(1,1000);
n1 = 50;
alpha = 5e-3;
gradW = W*(H*H') - A*H';
gradH = (W'*W)*H - W'*A;
Wn = max(W-alpha*gradW,0);
Hn = max(H-alpha*gradH,0);
newobj = sum([gradW; gradH'].^2,'all');
obj = sum((A-Wn*Hn).^2,'all');
% newobj = norm([gradW; gradH'],'fro');
% obj = 0.5*(norm(A-W*H,'fro')^2);


% main loop
% while newobj-obj > 0.1*(sum(sum(gradW.*(Wn-W)))+sum(sum(gradH.*(Hn-H))))
%     gradW = W*(H*H') - A*H';
%     gradH = (W'*W)*H - W'*A;
%     Wn = max(W-alpha*gradW,0);
%     Hn = max(H-alpha*gradH,0);
%     newobj = sum((A-W*H).^2,'all');
% %     newobj = 0.5*(norm(A-Wn*Hn,'fro')^2);
% end

for i = 1 : n1
    gradW = Wn*(Hn*Hn') - A*Hn';
    gradH = (Wn'*Wn)*Hn - Wn'*A;
    Wn = max(Wn-alpha*gradW,0);
    Hn = max(Hn-alpha*gradH,0);
    newobj = sum((A-Wn*Hn).^2,'all');
    Fnorm(1,i) = newobj;
%     newobj = 0.5*(norm(A-Wn*Hn,'fro')^2);
end

for i = n1+1 : 1000
    W = Wn.*(A*Hn')./(Wn*(Hn*Hn'));
    H = (Wn'*A).*Hn./((Wn'*Wn)*Hn);
    Wn = W;
    Hn = H;
    newobj = sum((A-Wn*Hn).^2,'all');
    Fnorm(1,i) = newobj;
%     newobj = 0.5*(norm(A-Wn*Hn,'fro')^2);
end
end
    
                
                
function [xiter,lm] = ASM(x,gfun,Hfun,A,b,W)
%% minimization using the active set method (Nocedal & Wright, Section 16.5)
% Solves f(x) --> min subject to Ax >= b
% x = initial guess, a column vector
TOL = 1e-10;
dim = length(x);
g = gfun(x);
H = Hfun(x);
iter = 0;
itermax = 20000;
m = size(A,1); % the number of constraints
% W = working set, the set of active constrains
I = (1:m)';
Wc = I; % the compliment of W
xiter = x;
while iter < itermax
    % compute step p: solve 0.5*p'*H*p + g'*p --> min subject to A(W,:)*p = 0
    AW = A(W,:); % LHS of active constraints
    % fix H if it is not positive definite
    ee = sort(eig(H),'ascend');
    if ee(1) < 1e-10
        lam = -ee(1) + 1;
    else
        lam = 0;
    end
    H = H + lam*eye(dim);
    if ~isempty(W)
        M = [H, -AW';AW,zeros(size(W,1))];
        RHS = [-g;zeros(size(W,1),1)];
    else
        M = H;
        RHS = -g;
    end
    aux = M\RHS;
    p = aux(1:dim);
    lm = aux(dim+1:end);
    if  norm(p) < TOL  % if step == 0
        if ~isempty(W)
            lm = AW'\g; % find Lagrange multipliers
            if min(lm) >= 0 % if Lagrange multipliers are positive, we are done
                % the minimizer is one of the corners
                fprintf('A local solution is found, iter = %d\n',iter);
                fprintf('x = [\n'); fprintf('%d\n',x);fprintf(']\n');
                break;
            else % remove the index of the most negative multiplier from W
                [lmin,imin] = min(lm);
                W = setdiff(W,W(imin));
                Wc = setdiff(I,W);
            end
        else
            fprintf('A local solution is found, iter = %d\n',iter);
            fprintf('x = [\n'); fprintf('%d\n',x);fprintf(']\n');
            break;
        end    
    else % if step is nonzero
        alp = 1;
        % check for blocking constraints
        Ap = A(Wc,:)*p;
        icand = find(Ap < -TOL);
        if ~isempty(icand)
            % find step lengths to all possible blocking constraints
            al = (b(Wc(icand)) - A(Wc(icand),:)*x)./Ap(icand);
            % find minimal step length that does not exceed 1
            [almin,kmin] = min(al);
            alp = min(1,almin);
        end
        x = x + alp*p;
        g = gfun(x);
        H = Hfun(x);
        if alp < 1
            W = [W;Wc(icand(kmin))];
            Wc = setdiff(I,W);
        end
    end
    iter = iter + 1;
    xiter = [xiter,x];
end
if iter == itermax
    fprintf('Stopped because the max number of iterations %d is performed\n',iter);
end
end

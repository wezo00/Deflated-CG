function[X,relres1,relres] = mult_rhs_DPCG(A, B, k, l, M)
% algorithm 5.3
% Inputs : A = nxn symmetric positive definite matrix
%          B = nxs matrix rhs 
%          k = nbr eigenvectors
%          l = nbr steps
%          M = LL', L = preconditioned matrix
%          
% Outputs : X, relres1 first system, relres last system

% solve first system with prec_defl_CG
[x, p, rel_res] = PCG(A,B(:,1),M,l);
X = x;
W = [];
relres1 = rel_res;
semilogy(relres1,"red");
hold on
% compute F^(1), G^(1)
Z = [W,p];
G = Z'*(A'*A)*Z;
F = Z'*A*Z;
% loop to solve systems 2 to the end
for s=2:width(B)
    [Y,~] = eigs(G,F,k,'smallestabs'); % generalized eigenvalue problem
    W = Z*Y;
    [x,P,rel_res] = prec_defl_CG(A,B(:,s),W,M,l); % solve next system with new W
    X = [X,x];
    relres = rel_res;
    semilogy(relres,"blue");
    % compute F^(s), G^(s)
    Z = [W,P];
    G = Z'*(A'*A)*Z;
    F = Z'*A*Z;
end
xlabel('number of iterations')
ylabel('true relative residual')
legend('PCG (s=1)','DeflatedCG (s=2:10)')
hold off

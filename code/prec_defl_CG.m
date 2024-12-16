function [x, p, rel_res] = prec_defl_CG(A, b, W, M, max_it)
% algorithm 3.6
% Goal : Approximate x such that Ax = b
% Inputs : A = nxn symmetric positive definite matrix
%          b = nx1 vector 
%          W = nxk matrix 
%          M = LL', L = preconditioned matrix
%          
% Outputs : x = [x_0, ..., x_j] = approximation of the solution of AX = b
%           p = [p_0, ..., p_j] = directions 
%           rel_res = [rel_res_0, ..., rel_res_j] = relative residual of
%                       each steps

n = length(b);
k = size(W, 2);

% Initialisation
x = zeros(n, max_it);
r = zeros(n, max_it);
z = zeros(n, max_it);
mu = zeros(k, max_it);
p = zeros(n, max_it);
D = zeros(1, max_it);
alpha = zeros(1, max_it);
beta = zeros(1, max_it);
rel_res = zeros(1, max_it);

% compute x_{-1} and r_{-1} :
x(:,1) = zeros(n,1);
r(:,1) = b - A*x(:,1);

% compute x0 and r0 with equation (3.12): 
x(:,1) = x(:,1) + W*((W'*A*W)\(W'*r(:,1)));
r(:,1) = b - A*x(:,1);

% solve W'AWmu0 = W'Az0 for mu, set p0 = -Wmu0 + z0
z(:,1) = M\r(:,1);
mu(:,1) = (W'*A*W)\(W'*A*z(:,1));
p(:,1) = -W*mu(:,1) + z(:,1);
D(1) = p(:,1)'*A*p(:,1);
rel_res(1) = norm(b - A*x(:,1), 2) / norm(b, 2);

j = 1;
while (rel_res(j) > 10^(-7))
    j = j + 1;
    alpha(j-1) = (r(:,j-1)' * z(:,j-1)) / (p(:,j-1)' * A * p(:,j-1));
    x(:,j) = x(:,j-1) + alpha(j-1) * p(:,j-1);
    r(:,j) = r(:,j-1) - alpha(j-1) * (A * p(:,j-1));
    z(:,j) = M \ r(:,j);
    beta(j-1) = (r(:,j)' * z(:,j)) / (r(:,j-1)' * z(:,j-1));
    mu(:,j) = (W'*A*W) \ (W'*A*z(:,j));
    p(:,j) = beta(j-1) * p(:,j-1) + z(:,j) - W * mu(:,j);
    D(j) = p(:,j)' * A * p(:,j);
    rel_res(j) = norm(b - A*x(:,j), 2) / norm(b, 2);
end

x = x(:,1:j);
p = p(:,1:j);
rel_res = rel_res(1:j);

end

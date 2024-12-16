function [x, p, rel_res] = PCG(A, b, M, l)
x_1 = zeros(size(b));

M_inv = inv(M);

x = x_1;
r = b - A*x;
p = M_inv * r;

D = p(:,1)'*A*p(:,1);
alpha = [];
beta = [];
rel_res = norm(b-A*x, 2)/norm(b, 2);
j = 1;

while (rel_res(j) > 10^(-7))
    j = j + 1;
    alpha = (r(:,j-1)' * M_inv*r(:,j-1))/(p(:,j-1)'*A*p(:,j-1));
    x = x + alpha*p(:,j-1);
    r = [r, r(:,j-1) - alpha*A*p(:,j-1)];
    beta = (r(:,j)'*M_inv*r(:,j))/(r(:,j-1)'*M_inv*r(:,j-1));
    p = [p, beta*p(:,j-1) + M_inv * r(:,j)];
    rel_res = [rel_res, norm(b-A*x, 2)/norm(b, 2)];
end
p = p(:,1:l);

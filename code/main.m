%% Example 1
close all
clear all
clc

N = 20; % test with N = 20, 40, 68
n = N^2;
A = gallery('poisson',N);
rng(5)
b = randn(n, 1);

% W = [w1, w2, w3, w4] are the eigenvectors corresponding to the four smallest
% eigenvalues of A
[~,largest_eigenval] = eigs(A,1,'largestabs');
[W,smallest_eigenval] = eigs(A,4,'smallestabs');

% CG 
[~, ~, rel_res] = PCG(A, b, eye(n), 0);
condition_number_A = largest_eigenval/smallest_eigenval(1,1);
% we compute the condition number like this because the function eigs of
% matlab gives us incorrect results

% Deflated-CG with W = [w1] and 
[~, ~, rel_res1] = prec_defl_CG(A, b, W(:,1), eye(n), 200);
condition_number_syst1 = largest_eigenval/smallest_eigenval(2,2);

% Deflated-CG with W = [w1, w2] 
[~, ~, rel_res2] = prec_defl_CG(A, b, W(:,1:2), eye(n), 200);
condition_number_syst2 = largest_eigenval/smallest_eigenval(3,3);

% Deflated-CG with W = [w1, w2, w3] 
[~, ~, rel_res3] = prec_defl_CG(A, b, W(:, 1:3), eye(n), 200);
condition_number_syst3 = largest_eigenval/smallest_eigenval(4,4);

% plot the results
figure (1)
semilogy(1:size(rel_res,2), rel_res, "red", 1:size(rel_res1,2), rel_res1, "magenta", 1:size(rel_res2,2), rel_res2, "blue", 1:size(rel_res3,2), rel_res3, "green")
xlabel('number of iterations')
ylabel('true relative residual')
legend('CG', 'Deflated-CG with W = [w1]', 'Deflated-CG with W = [w1, w2]', 'Deflated-CG with W = [w1, w2, w3]')


%% Example 3 (instead of example 2)
close all
clear all
clc

% load data
[A, rows, cols, entries, rep, field, symm] = mm_to_msm("1138_bus.mtx");
L = ichol(A);
s = 10;
k = 5;
l = 20;
rng(7)
B = randn(rows,s);

% run and plot
[X,relres1,relres_final] = mult_rhs_DPCG(A,B,k,l,L*L');

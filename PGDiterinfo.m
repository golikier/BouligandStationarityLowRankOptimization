function [s, U, V, f, time] = PGDiterinfo(r, s0, U0, V0, f0, f1, a0, b, c, iter)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements monotone PGD [OW25, Algorithm 4.2 with l = 0 or p = 1]
% on the real determinantal variety.
% Input:
%   - a positive integer r;
%   - a row vector s0 of at most r positive real numbers in decreasing order;
%   - an m-by-length(s0) matrix U0 having orthonormal columns, where m > r;
%   - an n-by-length(s0) matrix V0 having orthonormal columns, where n > r;
%   - functions f0 and f1 that, given an m-by-n matrix, return respectively
%     the objective function and its gradient at that matrix;
%   - a positive real number a0;
%   - b and c in (0, 1);
%   - a positive integer iter, the number of iterations to be performed.
% Output:
%   - cells s, U, and V containing all iterates;
%   - a vector f containing the function value at each iteration;
%   - a vector time containing the running time at each iteration.
iter = iter+1;
s = cell([iter 1]);
U = cell([iter 1]);
V = cell([iter 1]);
f = zeros(iter, 1);
time = zeros(iter, 1);
s{1} = s0;
U{1} = U0;
V{1} = V0;
tic
X = (U0.*s0)*V0';
f(1) = f0(X);
for i = 2:iter
    G = f1(X);
    a = a0;
    [U{i}, s{i}, V{i}] = svds(X-a*G, r, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300); % large-scale (truncated) SVD
    s{i} = diag(s{i})';
    s{i} = s{i}(s{i} > 0);
    r_new = length(s{i});
    U{i} = U{i}(:, 1:r_new);
    V{i} = V{i}(:, 1:r_new);
    X_new = (U{i}.*s{i})*V{i}';
    f(i) = f0(X_new);
    while f(i) > f(i-1) + c*trace(G'*(X_new-X))
        a = a*b;
        [U{i}, s{i}, V{i}] = svds(X-a*G, r, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300); % large-scale (truncated) SVD
        s{i} = diag(s{i})';
        s{i} = s{i}(s{i} > 0);
        r_new = length(s{i});
        U{i} = U{i}(:, 1:r_new);
        V{i} = V{i}(:, 1:r_new);
        X_new = (U{i}.*s{i})*V{i}';
        f(i) = f0(X_new);
    end
    X = X_new;
    time(i) = toc;
end
end
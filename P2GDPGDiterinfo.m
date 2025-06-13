function [s, U, V, f, time, R] = P2GDPGDiterinfo(r, s0, U0, V0, f0, f1, g0, g1, a, b, c, Delta, iter)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements P2GD-PGD [OGA24, Definition 7.1].
% Input:
%   - a positive integer r;
%   - a row vector s0 of at most r positive real numbers in decreasing order;
%   - an m-by-length(s0) matrix U0 having orthonormal columns, where m > r;
%   - an n-by-length(s0) matrix V0 having orthonormal columns, where n > r;
%   - functions f0 and f1 that, given an m-by-n matrix, return respectively
%     the objective function and its gradient at that matrix;
%   - functions g0 and g1 that, given (L, R) with L m-by-k, R n-by-k, and
%     k <= r, return respectively the objective function and its gradient
%     at L*R';
%   - a positive real number a;
%   - b and c in (0, 1);
%   - a positive real number Delta;
%   - a positive integer iter, the number of iterations to be performed.
% Output:
%   - cells s, U, and V containing all iterates;
%   - a vector f containing the function value at each iteration;
%   - a vector time containing the running time at each iteration;
%   - the number R of PGD iterations.
iter = iter+1;
s = cell([iter 1]);
U = cell([iter 1]);
V = cell([iter 1]);
f = zeros(iter, 1);
time = zeros(iter, 1);
R = 0;
s{1} = s0;
U{1} = U0;
V{1} = V0;
f(1) = g0(U0.*s0, V0);
tic
for i = 2:iter
    r_now = length(s{i-1});
    if s{i-1}(r_now) > Delta
        [s{i}, U{i}, V{i}, f(i), ~] = P2GDmap(r, r_now, s{i-1}, U{i-1}, V{i-1}, g0, g1, a, b, c);
    else
        R = R+1;
        [s{i}, U{i}, V{i}, f(i)] = PGDmap(r, s{i-1}, U{i-1}, V{i-1}, f0, f1, a, b, c);
    end
    time(i) = toc;
end
end
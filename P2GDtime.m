function [s, U, V, f, time] = P2GDtime(r, s, U, V, g0, g1, a, b, c, f_tol)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements P2GD [SU15, Algorithm 3].
% Input:
%   - a positive integer r;
%   - a row vector s of at most r positive real numbers in decreasing order;
%   - an m-by-length(s) matrix U having orthonormal columns, where m > r;
%   - an n-by-length(s) matrix V having orthonormal columns, where n > r;
%   - functions g0 and g1 that, given (L, R) with L m-by-k, R n-by-k, and
%     k <= r, return respectively the objective function and its gradient
%     at L*R';
%   - a positive real number a;
%   - b and c in (0, 1);
%   - a positive real number f_tol.
% Output:
%   - the first (s, U, V) such that g0(U.*s, V) <= f_tol;
%   - the value f of g0 at (U.*s, V);
%   - the running time required by P2GD to generate (s, U, V).
% The method is stopped if its running time exceeds 60 seconds.
% This function computes only what is necessary to generate (s, U, V).
tic
f = g0(U.*s, V);
while f > f_tol && toc <= 60
    [s, U, V, f, ~] = P2GDmap(r, length(s), s, U, V, g0, g1, a, b, c);
end
time = toc;
end
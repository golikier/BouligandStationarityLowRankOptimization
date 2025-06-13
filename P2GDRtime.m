function [s, U, V, f, time] = P2GDRtime(r, s, U, V, g0, g1, a, b, c, Delta, f_tol)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements P2GDR [OGA24, Definition 6.1].
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
%   - a positive real number Delta;
%   - a positive real number f_tol.
% Output:
%   - the first (s, U, V) such that g0(U.*s, V) <= f_tol;
%   - the value f of g0 at (U.*s, V);
%   - the running time required by P2GDR to generate (s, U, V).
% This function computes only what is necessary to generate (s, U, V).
tic
f = g0(U.*s, V);
while f > f_tol
    [s, U, V, f] = P2GDRmap(r, s, U, V, g0, g1, a, b, c, Delta);
end
time = toc;
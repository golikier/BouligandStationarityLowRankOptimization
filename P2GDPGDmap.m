function [s, U, V, f] = P2GDPGDmap(r, s, U, V, f0, f1, g0, g1, a, b, c, Delta)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements the P2GD-PGD map [OGA24, Algorithm 7.2].
% Input:
%   - a positive integer r;
%   - a row vector s of at most r positive real numbers in decreasing order;
%   - an m-by-length(s) matrix U having orthonormal columns, where m > r;
%   - an n-by-length(s) matrix V having orthonormal columns, where n > r;
%   - functions f0 and f1 that, given an m-by-n matrix, return respectively
%     the objective function and its gradient at that matrix;
%   - functions g0 and g1 that, given (L, R) with L m-by-k, R n-by-k, and
%     k <= r, return respectively the objective function and its gradient
%     at L*R';
%   - a positive real number a;
%   - b and c in (0, 1);
%   - a positive real number Delta.
% Output:
%   - the vector s and matrices U and V obtained by the P2GD-PGD map;
%   - the value f of f0 at (U.*s)*V'.
r_now = length(s);
if s(r_now) > Delta
    [s, U, V, f, ~] = P2GDmap(r, r_now, s, U, V, g0, g1, a, b, c);
else
    [s, U, V, f] = PGDmap(r, s, U, V, f0, f1, a, b, c);
end
end
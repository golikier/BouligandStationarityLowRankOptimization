function [s, U, V, f] = PGDmap(r, s, U, V, f0, f1, a, b, c, mu)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements the PGD map [OW25, Algorithm 4.1], which
% performs a backtracking projected line search along the negative
% gradient, on the real determinantal variety.
% Input:
%   - a positive integer r;
%   - a row vector s of at most r positive real numbers in decreasing order;
%   - an m-by-length(s) matrix U having orthonormal columns, where m > r;
%   - an n-by-length(s) matrix V having orthonormal columns, where n > r;
%   - functions f0 and f1 that, given an m-by-n matrix, return respectively
%     the objective function and its gradient at that matrix;
%   - a positive real number a;
%   - b and c in (0, 1);
%   - an optional real number mu greater than f0((U.*s)*V'), whose default
%     value is f0((U.*s)*V').
% Output:
%   - the vector s and matrices U and V obtained by the PGD map;
%   - the value f of f0 at (U.*s)*V'.
X = (U.*s)*V';
G = f1(X);
if nargin == 9
    mu = f0(X);
end
[U, s, V] = svds(X-a*G, r, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300); % large-scale (truncated) SVD
s = diag(s)';
s = s(s > 0);
r_new = length(s);
U = U(:, 1:r_new);
V = V(:, 1:r_new);
Y = (U.*s)*V';
f = f0(Y);
while f > mu + c*trace(G'*(Y-X))
    a = a*b;
    [U, s, V] = svds(X-a*G, r, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300); % large-scale (truncated) SVD
    s = diag(s)';
    s = s(s > 0);
    r_new = length(s);
    U = U(:, 1:r_new);
    V = V(:, 1:r_new);
    Y = (U.*s)*V';
    f = f0(Y);
end
end
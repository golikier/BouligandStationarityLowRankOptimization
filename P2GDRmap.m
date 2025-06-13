function [s_new, U_new, V_new, f_new] = P2GDRmap(r, s, U, V, g0, g1, a, b, c, Delta)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements the P2GDR map [OGA24, Algorithm 6.1], following
% closely the pseudocode in [OA23, Algorithm 6].
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
%   - a positive real number Delta.
% Output:
%   - the vector s_new and matrices U_new and V_new obtained by the P2GDR map;
%   - the value f_new of g0 at (U_new.*s_new, V_new).
r_now = length(s);
[s_new, U_new, V_new, f_new, ~] = P2GDmap(r, r_now, s, U, V, g0, g1, a, b, c);
r_diff = r_now-sum(s > Delta); % difference between the rank and the Delta-rank
if r_diff > 0 % rank reduction mechanism
    r_R = r_now;
    for j = 1:r_diff
        r_R = r_R-1;
        [s_R, U_R, V_R, f_R, ~] = P2GDmap(r, r_R, s(1:r_R), U(:, 1:r_R), V(:, 1:r_R), g0, g1, a, b, c);
        if f_R < f_new
            s_new = s_R;
            U_new = U_R;
            V_new = V_R;
            f_new = f_R;
        end
    end
end
end
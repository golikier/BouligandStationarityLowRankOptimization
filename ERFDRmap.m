function [s_new, U_new, V_new, f_new] = ERFDRmap(m, n, r, s, U, V, g0, g1, a, b, c, subclass, Delta)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements four subclasses of the ERFDR map [OA24, Algorithm 4.2]:
%   - the RFDR map [OA23, Algorithm 2] if subclass = 0;
%   - the CRFDR map with the ith cone from [OA24, Table 6.1] if subclass = i.
% The implementation follows closely the pseudocode in [OA24, Algorithm 7.2].
% Input:
%   - positive integers m, n, and r such that r < min{m, n};
%   - a row vector s of at most r positive real numbers in decreasing order;
%   - an m-by-length(s) matrix U having orthonormal columns;
%   - an n-by-length(s) matrix V having orthonormal columns;
%   - functions g0 and g1 that, given (L, R) with L m-by-k, R n-by-k, and
%     k <= r, return respectively the objective function and its gradient
%     at L*R';
%   - a positive real number a;
%   - b and c in (0, 1);
%   - subclass in {0, 1, 2, 3} for choosing between RFDR and CRFDR with
%     each of the three cones from [OA24, Table 6.1];
%   - a positive real number Delta.
% Output:
%   - the vector s and matrices U and V obtained by the ERFDR map;
%   - the value f_new of g0 at (U.*s, V).
r_now = length(s); % rank of the input
[s_new, U_new, V_new, f_new] = ERFDmap(m, n, r, r_now, s, U, V, g0, g1, a, b, c, subclass);
if r_now == r && s(r) <= Delta % rank reduction mechanism
    [s_R, U_R, V_R, f_R] = ERFDmap(m, n, r, r-1, s(1:r-1), U(:, 1:r-1), V(:, 1:r-1), g0, g1, a, b, c, subclass);
    if f_R < f_new
        s_new = s_R;
        U_new = U_R;
        V_new = V_R;
        f_new = f_R;
    end
end
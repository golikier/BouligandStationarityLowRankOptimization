function [s, U, V, f, B] = P2GDmap(r, r_now, s, U, V, g0, g1, a, b, c)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements the P2GD map [OGA24, Algorithm 5.1], following
% closely the pseudocode in [OA23, Algorithm 5].
% Input:
%   - a positive integer r;
%   - a nonnegative integer r_now <= r;
%   - a length-r_now row vector s of positive real numbers in decreasing order;
%   - an m-by-r_now matrix U having orthonormal columns, where m > r;
%   - an n-by-r_now matrix V having orthonormal columns, where n > r;
%   - functions g0 and g1 that, given (L, R) with L m-by-k, R n-by-k, and
%     k <= r, return respectively the objective function and its gradient
%     at L*R';
%   - a positive real number a;
%   - b and c in (0, 1).
% Output:
%   - the vector s and matrices U and V obtained by the P2GD map;
%   - the value f of g0 at (U.*s, V);
%   - a measure of B-stationarity, B, before applying the P2GD map.
%% Projection of the negative gradient onto the tangent space
G = -g1(U.*s, V);
G_U = U'*G;
G_V = G*V;
G_UV = G_U*V;
G_UUV = U*G_UV;
[U_perp, R_U, r_U] = qrcp(G_V-G_UUV, r_now);
[V_perp, R_V, r_V] = qrcp(G_U'-V*G_UV', r_now);
% Squared norm of the projection of the negative gradient onto the tangent space
B = norm(G_UV, 'fro')^2+norm(R_U, 'fro')^2+norm(R_V, 'fro')^2;
%% Backtracking projected line search
f_now = g0(U.*s, V);
if r_now == r % the tangent space is the tangent cone
    [s, U, V, f] = Backtrack(s, [U U_perp], [V V_perp], R_V', R_U, zeros(r_U, r_V), a);
else % the tangent space is a proper subset of the tangent cone
    G_perp = G-U*G_U+(G_UUV-G_V)*V';
    if norm(G_perp, 'fro') == 0
        [s, U, V, f] = Backtrack(s, [U U_perp], [V V_perp], R_V', R_U, zeros(r_U, r_V), a);
    else
        [U_perp_, s_perp, V_perp_] = svds(G_perp, r-r_now, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300); % large-scale (truncated) SVD
        s_perp = diag(s_perp)';
        s_perp = s_perp(s_perp > 0);
        r_perp = length(s_perp);
        U_perp_ = U_perp_(:, 1:r_perp);
        V_perp_ = V_perp_(:, 1:r_perp);
        B = B+s_perp*s_perp';
        R_U_1 = U_perp'*U_perp_;
        [U_perp_bar, R_U_2, r_U_bar] = qrcp(U_perp_-U_perp*R_U_1, r_perp);
        R_V_1 = V_perp'*V_perp_;
        [V_perp_bar, R_V_2, r_V_bar] = qrcp(V_perp_-V_perp*R_V_1, r_perp);
        [s, U, V, f] = Backtrack(s, [U U_perp U_perp_bar], [V V_perp V_perp_bar], [R_V' zeros(r_now, r_V_bar)], [R_U; zeros(r_U_bar, r_now)], ([R_U_1; R_U_2].*s_perp)*[R_V_1; R_V_2]', a);
    end
end
    function [Q, R, r_] = qrcp(A, r_)
        [Q, R, p] = qr(A, 'econ', 'vector');
        p(p) = 1:r_;
        r_ = sum(abs(diag(R)) > 1e-14);
        Q = Q(:, 1:r_);
        R = R(1:r_, :);
        R = R(:, p);
    end
    function [s, U, V, f] = Backtrack(s, U, V, R_12, R_21, R_22, a)
        [U_hat, s_hat, V_hat] = svds([diag(s)+a*G_UV a*R_12; a*R_21 a*R_22], r, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300);
        s_hat = diag(s_hat)';
        s_hat = s_hat(s_hat > 0);
        r_hat = length(s_hat);
        U_hat = U_hat(:, 1:r_hat);
        V_hat = V_hat(:, 1:r_hat);
        f = g0(U*(U_hat.*s_hat), V*V_hat);
        while f > f_now-c*a*B
            a = a*b;
            [U_hat, s_hat, V_hat] = svds([diag(s)+a*G_UV a*R_12; a*R_21 a*R_22], r, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300);
            s_hat = diag(s_hat)';
            s_hat = s_hat(s_hat > 0);
            r_hat = length(s_hat);
            U_hat = U_hat(:, 1:r_hat);
            V_hat = V_hat(:, 1:r_hat);
            f = g0(U*(U_hat.*s_hat), V*V_hat);
        end
        s = s_hat;
        U = U*U_hat;
        V = V*V_hat;
    end
B = sqrt(B);
end
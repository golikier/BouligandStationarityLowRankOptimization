function [s, U, V, f] = ERFDmap(m, n, r, r_now, s, U, V, g0, g1, a, b, c, subclass)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements four subclasses of the ERFD map [OA24, Algorithm 3.1]:
%   - the RFD map if subclass = 0;
%   - the CRFD map with the ith cone from [OA24, Table 6.1] if subclass = i.
% The implementation follows closely the pseudocode in [OA24, Algorithm 7.1].
% Input:
%   - positive integers m, n, and r such that r < min{m, n};
%   - a nonnegative integer r_now <= r;
%   - a length-r_now row vector s of positive real numbers in decreasing order;
%   - an m-by-r_now matrix U having orthonormal columns;
%   - an n-by-r_now matrix V having orthonormal columns;
%   - functions g0 and g1 that, given (L, R) with L m-by-k, R n-by-k, and
%     k <= r, return respectively the objective function and its gradient
%     at L*R';
%   - a positive real number a;
%   - b and c in (0, 1);
%   - subclass in {0, 1, 2, 3} for choosing between RFD and CRFD with
%     each of the three cones from [OA24, Table 6.1].
% Output:
%   - the vector s and matrices U and V obtained by the ERFD map;
%   - the value f of g0 at (U.*s, V).
f_now = g0(U.*s, V);
G = -g1(U.*s, V);
if r_now == r
    G1 = G'*U;
    G2 = G*V;
    B1 = norm(G1, 'fro');
    B2 = norm(G2, 'fro');
    if B1 >= B2
        X2 = V.*s;
        B = B1*B1;
        V = X2+a*G1;
        f = g0(U, V);
        while f > f_now - c*a*B
            a = a*b;
            V = X2+a*G1;
            f = g0(U, V);
        end
        [V, s, U_hat] = svd(V, 'econ');
        s = diag(s)';
        s = s(s > 0);
        r_now = length(s);
        U_hat = U_hat(:, 1:r_now);
        V = V(:, 1:r_now);
        U = U*U_hat;
    else
        X1 = U.*s;
        B = B2*B2;
        U = X1+a*G2;
        f = g0(U, V);
        while f > f_now - c*a*B
            a = a*b;
            U = X1+a*G2;
            f = g0(U, V);
        end
        [U, s, V_hat] = svd(U, 'econ');
        s = diag(s)';
        s = s(s > 0);
        r_now = length(s);
        U = U(:, 1:r_now);
        V_hat = V_hat(:, 1:r_now);
        V = V*V_hat;
    end
else
    if subclass == 0
        G1 = G'*U;
        G2 = G*V;
        B1 = norm(G1, 'fro');
        B2 = norm(G2, 'fro');
        [U_bar, s_bar, V_bar] = svds(G-U*G1'-(G2-U*(G1'*V))*V', r-r_now, 'largest', 'Tolerance', 1e-14, 'MaxIterations', 300); % large-scale (truncated) SVD
        s_bar = diag(s_bar)';
        s_bar = s_bar(s_bar > 0);
        r_bar = length(s_bar);
        U_bar = U_bar(:, 1:r_bar);
        V_bar = V_bar(:, 1:r_bar);
        if B1 >= B2
            B = B1*B1 + s_bar*s_bar';
            X2 = V.*s;
            V = X2+a*G1;
            f = g0([U U_bar], [V a*V_bar.*s_bar]);
            while f > f_now - c*a*B
                a = a*b;
                V = X2+a*G1;
                f = g0([U U_bar], [V a*V_bar.*s_bar]);
            end
            [V, s, U_hat] = svd([V a*V_bar.*s_bar], 'econ');
            s = diag(s)';
            s = s(s > 0);
            r_hat = length(s);
            U_hat = U_hat(:, 1:r_hat);
            V = V(:, 1:r_hat);
            U = [U U_bar]*U_hat;
        else
            B = B2*B2 + s_bar*s_bar';
            X1 = U.*s;
            U = X1+a*G2;
            f = g0([U a*U_bar.*s_bar], [V V_bar]);
            while f > f_now - c*a*B
                a = a*b;
                U = X1+a*G2;
                f = g0([U a*U_bar.*s_bar], [V V_bar]);
            end
            [U, s, V_hat] = svd([U a*U_bar.*s_bar], 'econ');
            s = diag(s)';
            s = s(s > 0);
            r_hat = length(s);
            U = U(:, 1:r_hat);
            V_hat = V_hat(:, 1:r_hat);
            V = [V V_bar]*V_hat;
        end
    else
        if subclass == 1
            [s_bar, imax] = max(abs(G));
            [~, jmax] = max(s_bar);
            imax = imax(jmax);
            s_bar = G(imax, jmax);
            U_bar = zeros(m, 1);
            U_bar(imax, 1) = 1;
            V_bar = zeros(n, 1);
            V_bar(jmax, 1) = 1;
        elseif subclass == 2
            [s_bar, imax] = max(vecnorm(G, 2, 2));
            U_bar = zeros(m, 1);
            U_bar(imax, 1) = 1;
            V_bar = G(imax, :)'/s_bar;
        elseif subclass == 3
            [s_bar, jmax] = max(vecnorm(G, 2, 1));
            U_bar = G(:, jmax)/s_bar;
            V_bar = zeros(n, 1);
            V_bar(jmax, 1) = 1;
        end
        B = s_bar*s_bar;
        [U, R1, P1] = qr([U U_bar], 'econ', 'vector');
        [V, R2, P2] = qr([V V_bar], 'econ', 'vector');
        Utilde = U*R1(:, P1);
        U_new = Utilde.*[s a*s_bar];
        Vtilde = V*R2(:, P2);
        f = g0(U_new, Vtilde);
        while f > f_now - c*a*B
            a = a*b;
            U_new = Utilde.*[s a*s_bar];
            f = g0(U_new, Vtilde);
        end
        [Utilde, s, Vtilde] = svd(R1(:, P1).*[s a*s_bar]*R2(:, P2), 'econ');
        s = diag(s)';
        s = s(s > 0);
        r_tilde = length(s);
        Utilde = Utilde(:, 1:r_tilde);
        Vtilde = Vtilde(:, 1:r_tilde);
        U = U*Utilde;
        V = V*Vtilde;
    end
end
end
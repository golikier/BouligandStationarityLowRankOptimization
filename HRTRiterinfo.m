function [L, R, f, NormGrad, lambda, time] = HRTRiterinfo(L0, R0, f0, f1, f2, gamma, gamma_c, eta, iter)
%% Description
% Author: Guillaume Olikier (2025-06-13)
% This function implements HRTR [LKB23, Algorithm 1] with the first lift
% phi of [LKB23, (1.1)], namely phi(L, R) := L*R', the hook from
% [LKB23, Example 3.11], and the Cauchy step.
% Input:
%   - L and R, hooked m-by-r and n-by-r real matrices, respectively;
%   - f0, a real-valued function of a matrix variable;
%   - f1, the gradient of f0;
%   - f2, a function that returns the directional derivative of f1 at its
%     first argument along its second argument;
%   - a positive real number gamma;
%   - gamma_c and eta in (0, 1);
%   - a positive integer iter, the number of iterations to be performed.
% Output:
%   - cells L and R containing all iterates;
%   - a vector f containing the function value at each iteration;
%   - a vector NormGrad containing the norm of the gradient at each iterate but the last;
%   - a vector lambda containing the smallest eigenvalue of the Hessian at each iterate but the last;
%   - a vector time containing the running time of each iteration.
%% Lifted cost function
g = @(L, R) f0(L*R');
g_1_1 = @(L, R) f1(L*R')*R;
g_1_2 = @(L, R) f1(L*R')'*L;
g_2_1 = @(L, R, dL, dR) f2(L*R', dL*R'+L*dR')*R + f1(L*R')*dR;
g_2_2 = @(L, R, dL, dR) f2(L*R', dL*R'+L*dR')'*L + f1(L*R')'*dL;
%% Initialization
NormGrad = zeros(iter, 1);
lambda = zeros(iter, 1);
iter = iter+1;
L = cell([iter 1]);
R = cell([iter 1]);
f = zeros(iter, 1);
time = zeros(iter, 1);
L{1} = L0;
R{1} = R0;
f(1) = g(L0, R0);
[m, r] = size(L0);
[n, ~] = size(R0);
dim1 = m*r;
dim2 = n*r;
dim = dim1+dim2;
tic
for i = 2:iter
    G_1_1 = g_1_1(L{i-1}, R{i-1});
    G_1_2 = g_1_2(L{i-1}, R{i-1});
    SquaredNormGrad = norm(G_1_1, 'fro')^2+norm(G_1_2, 'fro')^2;
    NormGrad(i-1) = sqrt(SquaredNormGrad);
    %% Smallest eigenvalue and associated eigenvector of the Hessian of the lifted cost function
    H = zeros(dim);
    for j = 1:m
        for k = 1:r
            dL_basis = zeros(m, r);
            dL_basis(j, k) = 1;
            H(:, (j-1)*r+k) = [reshape(g_2_1(L{i-1}, R{i-1}, dL_basis, zeros(n, r))', [dim1, 1]) ; reshape(g_2_2(L{i-1}, R{i-1}, dL_basis, zeros(n, r))', [dim2, 1])];
        end
    end
    for j = 1:n
        for k = 1:r
            dR_basis = zeros(n, r);
            dR_basis(j, k) = 1;
            H(:, dim1+(j-1)*r+k) = [reshape(g_2_1(L{i-1}, R{i-1}, zeros(m, r), dR_basis)', [dim1, 1]) ; reshape(g_2_2(L{i-1}, R{i-1}, zeros(m, r), dR_basis)', [dim2, 1])];
        end
    end
    [e, lambda(i-1)] = eigs(H, 1, 'smallestreal', 'Tolerance', 1e-14, 'MaxIterations', 30);
    e = e/norm(e);
    E_1 = reshape(e(1:dim1), [r m])';
    E_2 = reshape(e(dim1+1:dim), [r n])';
    %% Trust-region step
    if trace(G_1_1'*E_1) + trace(G_1_2'*E_2) > 0
        E_1 = -E_1;
        E_2 = -E_2;
    end
    if lambda(i-1) < 0 && SquaredNormGrad < -lambda(i-1)^3
        delta = -gamma*lambda(i-1);
        U_1 = E_1;
        U_2 = E_2;
    else
        delta = gamma*NormGrad(i-1);
        U_1 = -G_1_1;
        U_2 = -G_1_2;
    end
    normU = sqrt(norm(U_1, 'fro')^2+norm(U_2, 'fro')^2);
    a = trace(g_2_1(L{i-1}, R{i-1}, U_1, U_2)'*U_1)+trace(g_2_2(L{i-1}, R{i-1}, U_1, U_2)'*U_2);
    b = trace(G_1_1'*U_1)+trace(G_1_2'*U_2);
    t = delta/normU;
    if a > 0
        t = min([t -b/a]);
    end
    L{i} = L{i-1}+t*U_1;
    R{i} = R{i-1}+t*U_2;
    decreaseModel = -b*t-0.5*t*t*a;
    if decreaseModel == 0
        rho = 1;
    else
        f(i) = g(L{i}, R{i});
        rho = (f(i-1)-f(i))/decreaseModel;
    end
    while rho < eta
        delta = gamma_c*delta;
        t = delta/normU;
        if a > 0
            t = min([t -b/a]);
        end
        L{i} = L{i-1}+t*U_1;
        R{i} = R{i-1}+t*U_2;
        decreaseModel = -b*t-0.5*t*t*a;
        if decreaseModel == 0
            rho = 1;
        else
            f(i) = g(L{i}, R{i});
            rho = (f(i-1)-f(i))/decreaseModel;
        end
    end
    %% Hook
    [U_L, R_L] = qr(L{i}, 'econ');
    [V, s, U] = svd(R{i}*R_L', 'econ');
    s = sqrt(diag(s)');
    U = U_L*U;
    L{i} = U.*s;
    R{i} = V.*s;
    time(i) = toc;
end
end
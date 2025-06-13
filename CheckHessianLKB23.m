%% Check the Hessian for the problem from [LKB23, section 2.2].
% Author: Guillaume Olikier (2024-06-22)
% Notation from [OGA24, section 8.2].
% Function Q
D = diag([1 0.5]);
DD = diag([1 0.25]);
Y_opt = diag([1 0]);
Q = @(Y) 0.5*norm(D*(Y-Y_opt), 'fro')^2;
%Q = @(Y) 0.5*trace(DD*(Y-Y_opt)*(Y'-Y_opt));
Q_1 = @(Y) DD*(Y-Y_opt);
% Function psi
psi = @(t) 0.25*(t*(t*(t*t-2)-4)-2);
psi_1 = @(t) t*(t*t-1)-1;
psi_2 = @(t) 3*t*t-1;
% Function f
f = @(X) Q(X(1:2,1:2)) + psi(X(3,3));
f_1 = @(X) blkdiag(Q_1(X(1:2,1:2)), psi_1(X(3,3)));
f_2 = @(X, dX) blkdiag(DD*dX(1:2,1:2), psi_2(X(3,3))*dX(3,3));
% Test
X = randn(3);
dX = randn(3);
CheckHessian(f, f_1, f_2, X, dX);
%% Check of the Hessian for the lifted problem of [LKB23, section 2.2].
% Lift [OGA24, Proposition A.1]
g = @(L, R) f(L*R');
g_1_1 = @(L, R) f_1(L*R')*R;
g_1_2 = @(L, R) f_1(L*R')'*L;
g_2_1 = @(L, R, dL, dR) f_2(L*R', dL*R'+L*dR')*R + f_1(L*R')*dR;
g_2_2 = @(L, R, dL, dR) f_2(L*R', dL*R'+L*dR')'*L + f_1(L*R')'*dL;
% Test lift
L = randn(3, 2);
R = randn(3, 2);
dL = randn(3, 2);
dR = randn(3, 2);
CheckHessianLift(g, g_1_1, g_1_2, g_2_1, g_2_2, L, R, dL, dR);
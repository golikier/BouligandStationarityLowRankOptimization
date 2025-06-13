function [t, r] = CheckHessian(f, f_1, f_2, X, dX)
%% Description
% Author: Guillaume Olikier (2024-06-22)
% Input:
%   - a real-valued function f of a matrix variable;
%   - the gradient f_1 of f;
%   - a function f_2 that returns the directional derivative of f_1 at its
%     first argument along its second argument;
%   - matrices X and dX.
% Output:
%   - a vector t of points between 1e-6 and 1;
%   - a vector r containing the difference between f(X+t*dX) and its
%     second-order Taylor expansion
%     f(X) + <f_1(X), dX> * t + 0.5 * <f_2(X, dX), dX> * t^2.
% This function also plots r as a function of t in log-log scale. This
% plot, in blue, must look like the green line which has slope 3.
%% Code
i = 500;
t_log = linspace(-6, 0, i);
t = exp(t_log*log(10));
r = zeros(1, i);
a0 = f(X);
a1 = trace(f_1(X)'*dX);
a2 = trace(f_2(X, dX)'*dX);
for j = 1:i
    r(j) = abs(f(X+t(j)*dX)-(a0+a1*t(j)+0.5*a2*t(j)*t(j)));
end
figure
box
grid
hold on
plot(t_log, 3*log(t)/log(10), 'g.-', 'MarkerSize', 6);
plot(t_log, log(r)/log(10), 'b.', 'MarkerSize', 6);
xlabel('log_{10}(t)');
ylabel('log_{10}(|f(X+tV)-(f(X)-\langle\nablaf(X),V\ranglet-0.5\langle\nabla^2f(X)(V),V\ranglet^2)|)');
set(gca, 'FontSize', 12);
end
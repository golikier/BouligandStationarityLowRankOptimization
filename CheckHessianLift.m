function [t, r] = CheckHessianLift(g, g_1_1, g_1_2, g_2_1, g_2_2, L, R, dL, dR)
%% Description
% Author: Guillaume Olikier (2024-06-22)
% Input:
%   - a real-valued function g of two matrix variables;
%   - the gradient g_1_1 of g with respect to its first variable;
%   - the gradient g_1_2 of g with respect to its second variable;
%   - a function g_2_1 that returns the directional derivative of g_1_1 at
%     its first two arguments along its last two arguments;
%   - a function g_2_2 that returns the directional derivative of g_1_2 at
%     its first two arguments along its last two arguments;
%   - matrices L, R, dL, and dR.
% Output:
%   - a vector t of points between 1e-6 and 1;
%   - a vector r containing the difference between g(L+t*dL, R+t*dR) and
%     its second-order Taylor expansion.
% This function also plots r as a function of t in log-log scale. This
% plot, in blue, must look like the green line which has slope 3.
%% Code
i = 500;
t_log = linspace(-6, 0, i);
t = exp(t_log*log(10));
r = zeros(1, i);
a0 = g(L, R);
a1 = trace(g_1_1(L, R)'*dL)+trace(g_1_2(L, R)'*dR);
a2 = trace(g_2_1(L, R, dL, dR)'*dL)+trace(g_2_2(L, R, dL, dR)'*dR);
for j = 1:i
    r(j) = abs(g(L+t(j)*dL, R+t(j)*dR)-(a0+a1*t(j)+0.5*a2*t(j)*t(j)));
end
figure
box
grid
hold on
plot(t_log, 3*log(t)/log(10), 'g.-', 'MarkerSize', 6);
plot(t_log, log(r)/log(10), 'b.', 'MarkerSize', 6);
xlabel('log_{10}(t)');
ylabel('log_{10}(|g((L,R)+t(dL,dR))-(g(L,R)-\langle\nablag(L,R),(dL,dR)\ranglet-0.5\langle\nabla^2g(L,R)(dL,dR),(dL,dR)\ranglet^2)|)');
set(gca, 'FontSize', 12);
end
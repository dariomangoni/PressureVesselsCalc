clear
close all
clc

% Define variables
E = 210e9;  % Replace with the value of E (D9)
t = 3e-3;  % Replace with the value of t (D13)
a = 1.2/2;  % Replace with the value of a (D18)
n = 1:200;  % Replace with the value of n (B11)
l = 0.8;  % Replace with the value of l (D7)
v = 0.3;  % Replace with the value of v (D10)

P_cr = (E * (t / a)) ./ (n.^2 - 1 + 0.5 * (pi * a / l)^2) .* (1 ./ (n.^2 .* (l / (pi * a))^2 + 1).^2 + (t^2 ./ (12 * a^2 * (1 - v)^2)) .* (n.^2 - 1 + (pi * a / l)^2).^2);
[P_cr_min, P_cr_min_i] = min(P_cr, [], 'all');
n_worst = n(P_cr_min_i);
disp(['Worst case with n: ', num2str(n_worst)]);
disp(['P_cr: ', num2str(P_cr_min/100000), ' bar']);

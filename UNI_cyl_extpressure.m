clear
close all
clc

%% Design by EN 13445-3 - Clause 8

% Rp02T: 0,2% proof strength at temperature T
% sigma_e: % elastic limit for shell
% S: safety factor
% circ_tol: circularity tolerance

MaterialS235JR;
S_oper = 1.5; % operating conditions (8.4.4-1) 
S_test = 1.1; % testing conditions (8.4.4-2)
ea = 4; % thickness in analysis
R = 600; % mm, mean radius
h = 192; % mm, height of the head
Lcyl = 400; % length of cylinder (tangent-to-tangent)
E = 210e3; % MPa, modulus of elasticity
nu = 0.3; % poisson


%% Calcs
PrPy_PmPy_cylinders = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0;
0, 0.125, 0.251, 0.375, 0.5, 0.605, 0.68, 0.72, 0.755, 0.78, 0.803, 0.822, 0.836, 0.849, 0.861, 0.87, 0.879, 0.887, 0.896, 0.905, 0.914, 0.917, 0.923, 0.929, 0.935, 0.941, 0.947, 0.953, 0.959];


sigma_e = mat.Rp02T; % (8.4.2-1)

L = Lcyl + 0.4*h; % (8.5.2-2) 
Py = sigma_e*ea/R; % (8.5.2-4) 

Z = pi*R/L; % eq (8.5.2-7) 
ncyl_range = 1:20;
eps = (1 ./ (ncyl_range.^2 - 1 + (Z^2)/2)) .* ...
      (1 ./ ((ncyl_range.^2 / Z^2) + 1).^2 + (ea^2 / (12 * R^2 * (1 - nu^2))) .* (ncyl_range.^2 - 1 + Z.^2).^2); % (8.5.2-6) 

Pm_range = E*ea*eps/R; % (8.5.2-5)

[Pm, ncyl_id] = min(Pm_range);
ncyl = ncyl_range(ncyl_id);

Pm_Py = Pm/Py;
assert(Pm_Py>=0, 'Out of plot range');
Pr_Py_fun = @(Pm_Py) interp1(PrPy_PmPy_cylinders(1,:), PrPy_PmPy_cylinders(2,:), Pm_Py, 'linear', PrPy_PmPy_cylinders(2,end));
Pr_Py = Pr_Py_fun(Pm_Py);

P_test = Pr_Py*Py/S_test;
P_oper = Pr_Py*Py/S_oper;

disp(['Max external pressure (operating): ', num2str(P_oper*10), ' bar'])
disp(['Max external pressure (testing): ', num2str(P_test*10), ' bar'])

%% Figures
figure('NumberTitle', 'off', 'Name', 'Pm'); hold on; grid on;
plot(ncyl_range, Pm_range);
plot(ncyl, Pm, 'ko');
xlabel('n_{cyl}'); ylabel('P_m');

figure('NumberTitle', 'off', 'Name', 'Pr_Py'); hold on; grid on;
Pm_Py_range = linspace(0, 8.0);
plot(Pm_Py_range, Pr_Py_fun(Pm_Py_range));
plot(Pm_Py, Pr_Py, 'ko')
xlabel('P_m/P_y'); ylabel('P_r/P_y');

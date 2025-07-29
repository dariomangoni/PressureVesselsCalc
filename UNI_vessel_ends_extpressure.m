clear
close all
clc

MaterialS235JR;

ea = 4;
Di = 1200;
S_oper = 1.5; % operating conditions (8.4.4-1) 
S_test = 1.1; % testing conditions (8.4.4-2)
sigma_e = Rp02T; % (8.4.2-1)
E = 210e3; % MPa, modulus of elasticity
R = Di; % as Tommasin

warning('Also vessel ends as intpressure must pass.')

if 2.4*sqrt(ea*R)>1.3*R
    error('Expand formula with info contained in 8.7.2')
end


PrPy_PmPy_heads = [0, 0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4, 4.5, 5.0, 5.5, 6, 6.5;
0, 0.09, 0.18, 0.255, 0.324, 0.386, 0.435, 0.479, 0.51, 0.533, 0.548, 0.565, 0.567, 0.57]; % Figure 8.5-5


Py = 2 * sigma_e * ea/R; %  (8.7.1-1)
Pm = 1.21*E*ea^2/R^2; %  (8.7.1-2)

Pm_Py = Pm/Py;
assert(Pm_Py>=0, 'Out of plot range');
Pr_Py_fun = @(Pm_Py) interp1(PrPy_PmPy_heads(1,:), PrPy_PmPy_heads(2,:), Pm_Py, 'linear', PrPy_PmPy_heads(2,end));
Pr_Py = Pr_Py_fun(Pm_Py);

P_test = Pr_Py*Py/S_test;
P_oper = Pr_Py*Py/S_oper;

disp(['Max external pressure (operating): ', num2str(P_oper*10), ' bar'])
disp(['Max external pressure (testing): ', num2str(P_test*10), ' bar'])
 

%% Figures

figure('NumberTitle', 'off', 'Name', 'Pr_Py'); hold on; grid on;
Pm_Py_range = linspace(0, 8.0);
plot(Pm_Py_range, Pr_Py_fun(Pm_Py_range));
plot(Pm_Py, Pr_Py, 'ko')
xlabel('P_m/P_y'); ylabel('P_r/P_y');
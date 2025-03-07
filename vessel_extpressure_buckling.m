clear
close all
clc

%% general parameters
t = 4e-3; % thickness
len = 0.6; % cylinder length
poiss = 0.3; % poisson
E = 190e9; % Young modulus
stress_yield = 200e6; % yield stress
ovalization = 0.5e-2;

%% thin-walled parameters
r = 600e-3; % radius (inner radius if p positive, outer radius if p negative)


%% Ross (generic)
p = 1e5; % pressure (relative), positive if internal

% Thin-walled circular cylinder
stress_hoop_thin_cylinder = p*r/t; % [Ross, pag. 15]
stress_long_thin_cylinder = p*r/2/t; % [Ross, pag. 15]

disp(['Hoop stress: ', num2str(stress_hoop_thin_cylinder*1e-6), ' MPa'])
disp(['Longitudinal stress: ', num2str(stress_long_thin_cylinder*1e-6), ' MPa'])

%% Von Mises buckling [Ross, eq. 3.2]
disp('#### Von Mises ###')
p_buckling_fun = @(eig_n) E*t/r./(eig_n.^2 - 1 + 0.5*(pi*r/len).^2) .* (1./(eig_n.^2*(len/pi/r).^2 + 1).^2 + t^2./(12.*r^2*(1-poiss).^2).*(eig_n.^2-1+(pi*r/len).^2).^2);

p_buckling_VM = min(p_buckling_fun(2:10));
% figure; plot(p_buckling_fun(1:40));

disp(['Buckling at pressure: ', num2str(p_buckling_VM*1e-5), ' bar'])


%% Windenburg and Trilling [Ross, pag 105]
disp('#### Windenburg and Trilling ###')
thinness_ratio = ((len/(2*r)).^2 / (t/(2*r)).^3).^(1/4) * sqrt(stress_yield/E);

if 1/thinness_ratio<0.4
    PKD = 1;
else
    PKD = (2.8 - 1.0)/(1.46 - 0.4) * (1/thinness_ratio - 0.4) + 1;
end
p_buckling_WT = p_buckling_VM/PKD;
disp(['Buckling with PKD at pressure: ', num2str(p_buckling_WT*1e-5), ' bar'])

figure('Name', 'Plastic Knockdown Factor (PKD)', 'NumberTitle', 'off'); hold on; grid on;
line([1.0 1.0 2.8], [0.0 0.4 1.46], 'Color', 'k');
xlim([0 3.0]); ylim([0 1.6]);
plot(PKD, 1/thinness_ratio, 'r*');

%% Annaratone
disp('#### Annaratone ###')
disp('Critical pressures: ')
thinness_factor = t/(2*r);
%%  infinite length
p_hoop_inf = 2*stress_yield*thinness_factor; % [Annaratone, eq. 4.52] same as [Ross, pag. 15]
p_buck_inf = 2*E/(1-poiss^2)*thinness_factor^3; % [Annaratone, eq. 4.48] Bresse equation

%%  ovalization
C = 6*ovalization/4/thinness_factor; % [Annaratone, eq. 4.69]

p_hoop_inf_oval = p_hoop_inf/(1 + C); % [Annaratone, ]
p_buck_inf_oval = p_buck_inf/(1 + C*p_buck_inf/p_hoop_inf); % [Annaratone, ]

%%  finite length
p_buck_fin_fun = @(eig_n) 2*E./((eig_n.^2-1).*(1+(2*eig_n*len/pi/2/r).^2 ).^2)*t/2/r + 2*E/3/(1-poiss^2).*(eig_n.^2-1-(2*eig_n.^2 - 1 - poiss)/(1 - (2*eig_n*len/pi/2/r).^2))*(t/2/r)^3; % [Annaratone, eq. 4.91]
p_buck_fin = min(p_buck_fin_fun(2:10));
% p_buckling_finite = p_buckling_finite_fun(2);

%% finite length + ovalization
C1 = 6*(1-sqrt(2*r*t)/len)*ovalization/4/thinness_factor; % [Annaratone, eq. 4.104]
C2 = max([0, 6*(1-(2*r)/len*sqrt(thinness_factor))*(1-0.08*(2*r)/len/sqrt(thinness_factor))/thinness_factor*ovalization/4]); % [Annaratone, eq. 4.109]

p_hoop_fin_oval_C1 = p_hoop_inf/(1 + C1); % [Annaratone, eq. 4.94]
p_buck_fin_oval_C1 = p_buck_fin/(1 + C1*p_buck_fin/p_hoop_inf); % [Annaratone, eq. 4.108 modified with C1]

p_hoop_fin_oval_C2 = p_hoop_inf/(1 + C2); % [Annaratone, eq. 4.94]
p_buck_fin_oval_C2 = p_buck_fin/(1 + C2*p_buck_fin/p_hoop_inf); % [Annaratone, eq. 4.108 modified with C1]

disp(' ------ Hoop ------')
disp([' inf | circ : ', num2str(p_hoop_inf*1e-5), ' bar'])
disp([' inf | oval : ', num2str(p_hoop_inf_oval*1e-5), ' bar'])
disp(['*fin | oval : ', num2str(p_hoop_fin_oval_C1*1e-5), ' bar'])
disp([' fin | oval : ', num2str(p_hoop_fin_oval_C2*1e-5), ' bar'])
disp(' ------ Buckling ------')
disp([' inf | circ : ', num2str(p_buck_inf*1e-5), ' bar'])
disp([' fin | circ : ', num2str(p_buck_fin*1e-5), ' bar'])
disp([' inf | oval : ', num2str(p_buck_inf_oval*1e-5), ' bar'])
disp([' fin | oval : ', num2str(p_buck_fin_oval_C1*1e-5), ' bar'])
disp(['*fin | oval : ', num2str(p_buck_fin_oval_C2*1e-5), ' bar'])
disp(' ------ Summary ------')
disp(['Hoop : ', num2str(p_hoop_inf*1e-5), ' bar'])
if(2*r/len*sqrt(1/thinness_factor)>=12.5)
    % do not consider ovalization
    p_buck_critical = p_buck_fin;
else
    % consider ovalization
    p_buck_critical = p_buck_fin_oval_C2;
end
disp(['Buck : ', num2str(p_buck_critical*1e-5), ' bar'])
disp([' ---- Max allowed pressure : ', num2str(min([p_buck_critical, p_hoop_inf])*1e-5), ' bar'])


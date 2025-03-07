clear
close all
clc

D = 1200e-3;
mu = 0.3; % poisson
len = 0.8;
E = 210e9;
n = 6;
t = 6e-3;

rho = 1/(n^2*(2*len/pi/D)^2 + 1);
lam1 = rho*(2- rho)/(1-rho)^2;
lam2 = rho*(3+mu+(1-mu^2)*rho);
lam3 = rho*(1+mu) - rho^2*(mu*(1+2*mu) + (1-mu^2)*(1-rho*mu)*(1 + (1+mu)/(1-mu)*rho));

% radial pressure only
p_12B = 1/3*(n^2-1 + (lam1*n^4 - lam2*n^2 + lam3)/(n^2-1))*2*E/(1-mu^2)*(t/D)^3 + 2*E*t/D/((n^2-1)*(n^2*(2*len/pi/D)^2 + 1)^2);
disp(['p_12B: ', num2str(p_12B*1e-5), ' bar'])

p_12D = 1/3*(n^2-1+(2*n^2-1-mu)/(n^2*(2*len/pi/D)^2-1))*2*E/(1-mu^2)*(t/D)^3 + 2*E*t/D/((n^2-1)*(n^2*(2*len/pi/D)^2+1)^2);
disp(['p_12D: ', num2str(p_12D*1e-5), ' bar'])

% radial + axial
mu1 = (2 * lam2)/2;
mu2 = 1 + lam3;
p_13 = (1/3*((n^2+(pi*D/2/len)^2)^2 -2*mu1*n^2 + mu2)*2*E/(1-mu^2)*(t/D)^3 + 2*E*(t/D)/(n^2*(2*len/pi/D)^2 + 1)^2)/(n^2-1+0.5*(pi*D/2/len)^2);
disp(['p_13: ', num2str(p_13*1e-5), ' bar'])

r = D/2;
poiss = mu;
p_buckling_fun = @(eig_n) E*t/r./(eig_n.^2 - 1 + 0.5*(pi*r/len).^2) .* (1./(eig_n.^2*(len/pi/r).^2 + 1).^2 + t^2./(12.*r^2*(1-poiss).^2).*(eig_n.^2-1+(pi*r/len).^2).^2);

disp(['p Ross: ', num2str(p_buckling_fun(n)*1e-5), ' bar'])


thinness_factor = t/(2*r);
p_critical_buckling_finite_fun = @(eig_n) 2*E./((eig_n.^2-1).*(1+(2*eig_n*len/pi/2/r).^2 ).^2)*t/2/r + 2*E/3/(1-poiss^2).*(eig_n.^2-1-(2*eig_n.^2 - 1 - poiss)/(1 - (2*eig_n*len/pi/2/r).^2))*(t/2/r)^3;

disp(['p Annaratone: ', num2str(p_critical_buckling_finite_fun(n)*1e-5), ' bar'])

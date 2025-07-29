clear
close all
clc

% Di: internal diameter of shell
% P: pressure
% ea: analysis thickness
% f: nominal design stress
% z: joint efficiency

MaterialS235JR;

Di = 1200; % mm
e = 4; % mm
f = fd; % MPa % TO VERIFY
z = 0.8;
P = 3/10;
E = 210e3; 
mass = 100; %kg
Fadd = 0;
L = 0.5;


assert(e/(Di+2*e)<=0.16, 'Sec 7.4.1 - Formula not valid for this sizes.')

ea = e;
Dm = Di + e;

P_allowed = 2*f*z*ea/Dm;

disp(['Allowable pressure: ', num2str(P_allowed*10), ' bar'])


% % 16.14.6 Cylinder under internal pressure (P > 0)
% D = Dm;
% F = mass*9.81 + Fadd + pi/4*D^2*P;
% M = 0; % TO DO
% sigma_max = (F*D + 4*M)/(pi*D^2*ea);
% sigma_min = (F*D - 4*M)/(pi*D^2*ea);
% if sigma_min < 0
%     sigma_c = -sigma_min;
% end
% sigma_p = P*Dm/2/ea;
% 
% assert(sigma_max<=f)
% 
% if sigma_min < 0
%     error('Case not covered. Go to 16.14.6 Cylinder under internal pressure (P > 0), point 6')
% end
% 
% % 16.14.8 Global longitudinal compressive stress limits
% sigma_e = Rp02T;
% if D/ea <= 0.06*E/sigma_e
%     sigma_c_all = f;
% else
%     omega = L/sqrt(0.5*D*ea);   % eq (16.14-16) 
%     if omega<1.7
%         Cx = 1.36 - 1.83/omega + 2.07/omega^2;
%     elseif omega<=0.25*D/ea
%         Cx = 1.0;
%     else
%         Cxb = 1.0; % see Table 16.14–1
%         Cx = max([1.0+0./Cxb*(1-4*omega*ea/D), 0.6]);
%     end
% 
%     K = 1.21*E*Cx*ea/(sigma_e*D);
%     Q = 16; % see Table 16.14–2 — Fabrication quality parameter Q 
%     deltaWk = sqrt(0.5*D*ea)/Q;
%     alpha = 0.62/(1+1.91*(deltaWk/ea)^1.44);
%     lambdax0 = 0.2;
%     lambdap = sqrt(2.5*alpha);
%     lambdax = sqrt(1/K);
% 
%     if lambdax<=lambdax0
%         Delta = 1;
%     elseif lambdax<lambdap
%         Delta = 1-0.6*(lambdax - lambdax0)/(lambdap - lambdax0);
%     else
%         Delta = alpha*K;
%     end
%     S = 1.5; % see Table 5.3.2.4-1,
%     sigma_c_all = Delta*psi*sigma_e/S;
% end


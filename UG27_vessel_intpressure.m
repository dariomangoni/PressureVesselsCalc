clear
close all
clc

%% from [Moss]

% P: internal pressure, psi
% Di, Do: inside/outside diameter, in.
% S: allowable or calculated stress, psi
% E: joint efﬁciency (0.35 - 1.0)
% L: crown radius, in.
% Ri, Ro: inside/outside radius, in.
% K, M: coefﬁcients (See Note 3)
% sx: longitudinal stress, psi
% sf: circumferential stress, psi
% Rm: mean radius of shell, in.
% t: thickness or thickness required of shell, head, or cone, in.
% r: knuckle radius, in.

% h: depth of head, in.
% RL: latitudinal radius of curvature, in.
% Rm: meridional radius of curvature, in.
% sf: latitudinal stress, psi
% sx: meridional stress, psi
% P: internal pressure, psi
% tc: head thickness, corroded, in 

unit_conversion;

%% Parameters
r = 0.1 * FROM_M_TO_INCHES;
t = 4e-3 * FROM_M_TO_INCHES;
Ri = 600e-3 * FROM_M_TO_INCHES;
E = 0.95;
S = 200e6 * FROM_PA_TO_PSI;

if t>0.5*Ri
    error('Unsupported case')
end


%% Shell
P_circ = S*E*t/(Ri+0.6*t);
P_long = 2*S*E*t/(Ri-0.4*t);
P = min([P_long, P_circ]);
disp(['[SHELL] Pressure limit: ', num2str(P/FROM_PA_TO_PSI*1e-5), ' bar']);

if P > 0.385*S*E
    error('Look for other formulas: ASME BPVC 1-2 (pag 413 PDF)')
end
if P > 1.25*S*E
    error('Look for other formulas: ASME BPVC UG27 (pag 80 PDF)')
end
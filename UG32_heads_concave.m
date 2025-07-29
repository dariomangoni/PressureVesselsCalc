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
r = 0.05 * FROM_M_TO_INCHES;
t = 4e-3 * FROM_M_TO_INCHES;
Ri = 600e-3 * FROM_M_TO_INCHES;
E = 0.95;
S = 200e6 * FROM_PA_TO_PSI;

if S > 70000
    error('See note 90 ASME BPVC')
end

% Head (tori)
L = 0.9 * 2 * Ri;
h = 192e-3 * FROM_M_TO_INCHES;
tc = t*0.95;

%% General
Di = 2*Ri;
K = 0.167*(2 + (Di/2/h)^2);
M = 0.25*(3 + sqrt(L/r));

%% Heads
R = Ri;
P_UG32 = S*E*tc/(.885*L + 0.1*tc); % ASME BPVC Sec VIII Div 1 - UG-32 eq 2 (criteria for tc/L>=0.002)
if tc/L<0.002
    M = 0.25*(3 + sqrt(L/r));
    P_14_3 = 2*S*E*tc/(L*M +0.2*tc); % ASME BPVC Sec VIII Div 1 - 1-4 eq 4
    P_14_4 = 2*S*E*tc/(M*L - tc*(M-0.2)); % ASME BPVC Sec VIII Div 1 - 1-4 eq 4
    disp(['[TORI] Pressure limit (BPVC 1-4 eq3): ', num2str(P_14_3/FROM_PA_TO_PSI*1e-5), ' bar']);
    disp(['[TORI] Pressure limit (BPVC 1-4 eq4): ', num2str(P_14_4/FROM_PA_TO_PSI*1e-5), ' bar']);
    if tc < 0.0005
        error('Unhandled case. See BPVC 1-4 (f)(1)')
    end
end
disp(['[TORI] Pressure limit (BPVC UG-32): ', num2str(P_UG32/FROM_PA_TO_PSI*1e-5), ' bar']);
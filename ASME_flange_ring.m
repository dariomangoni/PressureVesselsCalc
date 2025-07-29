clear
close all
clc

unit_conversion;

% A: flange O.D., in.
% Ab: cross-sectional area of bolts, in.2
% Am: total required cross-sectional area of bolts, in.2
% a: nominal bolt diameter, in.

% B1: flange I.D., in.
% Bs: bolt spacing, in.
% b: effective gasket width, in.
% bo: gasket seating width, in.
% C: bolt circle diameter, in.
% d: hub shape factor
% d1: bolt hole diameter, in.
% E, hD, hG, hT, R: radial distances, in.
% e: hub shape factor
% F: hub shape factor for integral-type flanges
% FL: hub shape factor for loose-type flanges
% 
% f: hub stress correction factor for integral flanges
% G: diameter at gasket load reaction, in.
% go: thickness of hub at small end, in.
% g1: thickness of hub at back of flange, in.
% H: hydrostatic end force, lb
% HD: hydrostatic end force on area inside of flange, lb
% HG: gasket load, operating, lb
% Hp: total joint-contact surface compression load, lb
% HT: pressure force on flange face, lb
% h: hub length, in.
% ho: hub factor
% MD: moment due to HD, in.-lb
% MG: moment due to HG, in.-lb
% Mo: total moment on flange, operating, in.-lb
% M0o: total moment on flange, seating
% MT: moment due to HT, in.-lb
% m: gasket factor
% mo: unit load, operating, lb
% mg: unit load, gasket seating, lb
% N: width of gasket, in.
% n: number of bolts
% v: Poisson’s ratio, 0.3 for steel
% P: design pressure, psi
% Ra: single bolt root area, in2 (from Table 3-3)
% Sa: allowable stress, bolt, at ambient temperature, psi
% Sb: allowable stress, bolt, at design temperature, psi
% Sfa: allowable stress, flange, at ambient temperature, psi
% Sfo: allowable stress, flange, at design temperature, psi
% SH: longitudinal hub stress, psi
% SR: radial stress in flange, psi
% ST: tangential stress in flange, psi
% T, U, Y, Z: K-factors (see Table 3-6)
% Tr, Ur, Yr: K-factors for reverse flanges
% t: flange thickness, in.
% tn: pipe wall thickness, in.
% V: hub shape factor for integral flanges
% VL: hub shape factor for loose flanges
% W: flange design bolt load, lb
% Wm1: required bolt load, operating, lb
% Wm2: required bolt load, gasket seating, lb
% w: width of raised face or gasket contact width, in. (See Table 3-5)
% y: gasket design seating stress, psi

% Ab: cross-sectional area of bolts, in.2
% Ag: actual joint-contact area of gasket, in.2
% b: effective gasket seating width, in.
% d: root diameter of threads, in.
% dm: pitch diameter of threads, in.
% G: diameter at location of gasket load reaction, in.
% M: external bending moment, in.-lb
% m: gasket factor
% N: gasket width, in.
% n: number of bolts
% Eb: modulus of elasticity of bolting material at temperature, psi
% Eg: modulus of elasticity of gasket material at temperature, psi
% P: internal pressure, psi
% Pe: equivalent pressure including external loads, psi
% Pr: radial load, lb
% PT: test pressure, psi
% F: restoring force of gasket (decreasing compression force) from initial bolting strain, lb
% Fbo: initial tightening force, lb
% Lb: effective length of bolt, mid nut to mid nut, in.
% W: total tightening force, lb
% Wml: H þ Hp: required bolt load, operating, lb
% Wm2: required bolt load, gasket seating, lb
% y: gasket unit seating load, psi
% H: total hydrostatic end force, lb
% HP: total joint-contact surface compression load, lb
% T: initial tightening torque required, ft-lb
% tg: thickness of gasket, in.
% tn: thickness of nut, in.
% K_seating: total friction factor between bolt/nut and nut/flange face
% w: width of ring joint gasket, in.

% from [Moss, Table 3.3, pag 158 PDF]
bolt_flange.Ra_table = [0.126, 0.202, 0.302, 0.419, 0.551, 0.693, 0.89, 1.054, 1.294, 1.515, 1.744, 2.049, 2.3, 3.02, 3.715, 4.618, 5.621];
bolt_flange.bolt_diam_name = {'1/2in', '5/8in', '3/4in', '7/8in', '1in', '1-1/8in', '1-1/4in', '1-3/8in', '1-1/2in', '1-5/8in', '1-3/4in', '1-7/8in', '2in', '2-1/4in', '2-1/2in', '2-3/4in', '3in'};
bolt_flange.bolt_diam_table = [1/2, 5/8, 3/4, 7/8, 1, 1+1/8, 1+1/4, 1+3/8, 1+1/2, 1+5/8, 1+3/4, 1+7/8, 2, 2+1/4, 2+1/2, 2+3/4, 3];
bolt_flange.R_table = [0.813 , 0.938 , 1.125 , 1.25 , 1.375 , 1.5 , 1.75 , 1.875 , 2 , 2.125 , 2.25 , 2.375 , 2.5 , 2.75 , 3.063 , 3.375 , 3.625];
bolt_flange.E_table = [0.625 , 0.75 , 0.813 , 0.938 , 1.063 , 1.125 , 1.25 , 1.375 , 1.5 , 1.625 , 1.75 , 1.875 , 2 , 2.25 , 2.375 , 2.625 , 2.875];
bolt_flange.Bs_min_table = [ 1.25 , 1.5 , 1.75 , 2.0625 , 2.25 , 2.5 , 2.8125 , 3.0625 , 3.25 , 3.5 , 3.75 , 4 , 4 , 4.5 , 5 , 5.5 , 6];

% from 
bolt_torquing.bolt_diam_table = [3/4, 7/8, 1, 1+1/8, 1+1/4, 1+3/8, 1+1/2, 1+5/8, 1+3/4, 1+7/8, 2, 2+1/4, 2+1/2, 2+3/4, 3, 3+1/4];
bolt_torquing.d_table = [0.6273, 0.7387, 0.8466, 0.9716, 1.0966, 1.2216, 1.3466, 1.4716, 1.5966, 1.7216, 1.8466, 2.0966, 2.3466, 2.5966, 2.8466, 3.0966];
bolt_torquing.dm_table = [0.6850, 0.8028, 0.9188, 1.0438, 1.1688, 1.2938, 1.4188, 1.5438, 1.6688, 1.7938, 1.9188, 2.1688, 2.4188, 2.6688, 2.9188, 3.1688];
bolt_torquing.tn_table = [0.7344, 0.8594, 0.9844, 1.1094, 1.2188, 1.3438, 1.4688, 1.5938, 1.7188, 1.8438, 1.9688, 2.2031, 2.4531, 2.7031, 2.9531, 3.1875];

%% Parameters
% Flange
D  = 1200e-3 * FROM_M_TO_INCHES; % inside diameter of cylinder, in.
G = 1250e-3 * FROM_M_TO_INCHES; % diameter at gasket load reaction, in.
P = 6e5 * FROM_PA_TO_PSI; % design pressure, psi
b = 10e-3 * FROM_M_TO_INCHES; % effective gasket width, in.
Sa = 640e6 * FROM_PA_TO_PSI; % allowable stress, bolt, at ambient temperature, psi
Sb = 640e6 * FROM_PA_TO_PSI; % allowable stress, bolt, at design temperature, psi
Sfo = 235e6 * FROM_PA_TO_PSI; % allowable stress, flange, at ambient temperature, psi
Sfa = 235e6 * FROM_PA_TO_PSI; % allowable stress, flange, at design temperature, psi
B = 1206e-3 * FROM_M_TO_INCHES; % flange I.D. (internal diameter), in.
t = 50e-3 * FROM_M_TO_INCHES; % flange thickness, in.
g1 = 10e-3 * FROM_M_TO_INCHES; % thickness of hub at back of flange, in.

% Gasket
v = 0.3; % Poisson's ratio, 0.3 for steel
K_seating = 0.15; % total friction factor between bolt/nut and nut/flange face, from [Moss, pag 181 PDF]
Eb = 200e6 * FROM_PA_TO_PSI; % modulus of elasticity of bolting material at temperature, psi
Eg = 10e3; % modulus of elasticity of gasket material at temperature, psi (pag 181 PDF)
tg = 6e-3 * FROM_M_TO_INCHES; % thickness of gasket, in.
Lb = 40e-3 * FROM_M_TO_INCHES; % effective length of bolt, mid nut to mid nut, in.
M = 0; % external moment
Pr = 0; % external force
assert(G>B, 'Gasket reactions diameter cannot be smaller than the flange inner diameter');

%% Ring Flange Design
% m, y: from [Moss, Table 3-4, pag 167 PDF]
Pe = 16*M/pi/G^3 + 4*Pr/pi/G^2 + P;
m = 0.25;
y = 0;
Hp = 2*b*pi*G*m*Pe;
H = G.^2*pi*Pe/4;
Wm1 = Hp + H;
Wm2 = b*pi*G*y;
Am = max([Wm2/Sa, Wm1/Sb]);


% Step 1 - Bolt size and number
n_tentative = round(D/4)*4; 
Ra_tentative = Am/n_tentative;
bolt_flange_id = find(bolt_flange.Ra_table>Ra_tentative, 1);
bolt_flange_id = 4;
assert(~isempty(bolt_flange_id), 'The bolt_flange table does not have a bolt so big.')
Ra = bolt_flange.Ra_table(bolt_flange_id);
n = round(ceil(Am/Ra)/4)*4; % always multiple of 4
n = 24; % OVERRIDE number of bolts if needed
a = bolt_flange.bolt_diam_table(bolt_flange_id);
disp(['Bolt num: ', num2str(n)]);
disp(['Bolt size: ', num2str(bolt_flange.bolt_diam_name{bolt_flange_id}), ' | ', num2str(a/FROM_M_TO_INCHES*1000), ' mm']);

% Step 2 - Bolt circle diameter
R = bolt_flange.R_table(bolt_flange_id);
E = bolt_flange.E_table(bolt_flange_id);
C = B + 2*g1 + 2*R;
A = C + 2*E;
disp(['Bolt circle diameter: ', num2str(C), ' in | ', num2str(C/FROM_M_TO_INCHES*1000), ' mm'])
disp(['Flange outer diameter: ', num2str(A), ' in | ', num2str(A/FROM_M_TO_INCHES*1000), ' mm'])
assert(A>G, 'Flange outer diameter cannot be smaller of gasket reaction circle');
assert(A>C, 'Flange outer diameter cannot be smaller of bold circle diameter');
assert(C>G, 'Bolts must not be inside the gasket circle diameter.');


% Step 3 - Minimum bolt spacing (not ASME requirement)
Bs = C/n;
Bs_min = bolt_flange.Bs_min_table(bolt_flange_id); % alternatively Bs_min (due to wrench clearance) = 2.0|2.5 * a
disp(['Bolt spacing: ', num2str(Bs), ' in | ', num2str(Bs/FROM_M_TO_INCHES*1000), ' mm'])
assert(Bs>Bs_min, 'Bolt spacing is below minimum.')

% Maximum bolt spacing
Bs_max = 2*a + t;
if Bs>Bs_max
    mo_mG_increase_factor = sqrt(Bs/Bs_max);
    disp('!!! Bolt spacing above maximum !!!!')
    disp(['  I''m increasing unit loads by a factor of ', num2str(mo_mG_increase_factor)]);
    disp('  but consider to increase the number of bolts instead.');
else
    mo_mG_increase_factor = 1.0;
end
% assert(Bs<Bs_max, 'Bolt spacing is above maximum bolt spacing.')

% Gasket width check
if y>0
    Nmin = Ab*Sa/(2*pi*y*G);
    assert(2*b>Nmin, 'Gasket width is below limits.') % VERIFY: b and N depends on flange type [Moss, Fig. 3-16, pag 180 PDF]
end

% Momentums
Wmt = Hp + H;
HD = pi*B^2*Pe/4;
HG = Wmt - H;
HT = H - HD;
hD = 0.5*(C-B);
hG = 0.5*(C-G);
hT = 0.5*(hD + hG);
MD = HD*hD;
MG = HG*hG;
MT = HT*hT;
Mo = MD + MG + MT;
disp(['Total momentum (operating): ', num2str(Mo), ' in-lb | ', num2str(Mo/FROM_M_TO_INCHES/FROM_N_TO_LBF), ' Nm'])
Mo = Mo * mo_mG_increase_factor;
disp(['Total momentum (operating, increased): ', num2str(Mo), ' in-lb | ', num2str(Mo/FROM_M_TO_INCHES/FROM_N_TO_LBF), ' Nm'])

Ab = min([Ra*n, n*pi*a^2/4]); % different pages report different way of calculating it
W = 0.5*(Am + Ab)*Sa; % can be more precise after bolt torque calculation
HG_seating = W;
hG_seating = hG;
MG_seating = HG_seating*hG_seating;
Mo_seating = MG_seating; % TO BE VERIFIED: I assume that during seating there is no pressure so no terms MD and MT
disp(['Total momentum (seating): ', num2str(Mo_seating), ' in-lb | ', num2str(Mo_seating/FROM_M_TO_INCHES/FROM_N_TO_LBF), ' Nm'])
Mo_seating = Mo_seating * mo_mG_increase_factor;
disp(['Total momentum (seating, increased): ', num2str(Mo_seating), ' in-lb | ', num2str(Mo_seating/FROM_M_TO_INCHES/FROM_N_TO_LBF), ' Nm'])

% Flange thickness
K = A/B;
U = K^2*(1+4.6052*(1+v)/(1-v)*log10(K)-1)/(1.0472*(K^2-1)*(K-1)*(1+v));
Y = (1-v^2)*U;
t_operating = sqrt(Mo*Y/Sfo/B);
t_seating = sqrt(Mo_seating*Y/Sfa/B);
disp(['Thickness (operating): ', num2str(t_operating), ' in | ', num2str(t_operating/FROM_M_TO_INCHES*1000), ' mm']);
disp(['Thickness (seating): ', num2str(t_seating), ' in | ', num2str(t_seating/FROM_M_TO_INCHES*1000), ' mm']);
t = max([t_operating, t_seating]);

% Bolt loads
bolt_torquing_id = find(bolt_torquing.bolt_diam_table>a, 1);
assert(~isempty(bolt_torquing_id), 'The bolt_torquing table does not have a bolt so big.')
d = bolt_torquing.d_table(bolt_torquing_id);
dm = bolt_torquing.dm_table(bolt_torquing_id);
Ab_seating = pi*d^2*n/4;
Ag = 2*pi*b*G;
DF = H/(1+(Ab_seating*Eb*tg)/(Ag*Eg*Lb));
Fbo = Hp + DF;
W = max([Fbo, Wm2]);
T = K_seating*W*dm/(12*n);
disp(['Bolt torque for seating: ', num2str(T), ' ft-lb | ', num2str(T/FROM_N_TO_LBF/FROM_M_TO_FEET), ' Nm'])


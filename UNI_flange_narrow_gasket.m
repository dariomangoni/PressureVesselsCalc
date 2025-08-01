clear
close all
clc

% Section 11.10 Full face flanges with metal to metal contact
% refers to 11.5.2 with some modifications

% A is the outside diameter of the flange or, where slotted holes extend to outside of flange, the diameter to bottom of slots; 
% AB is the total cross-sectional area of bolts at the section of least bolt diameter; 
% AB,min is the total required cross-sectional area of bolts; 
% A2 is the outside diameter of the contact face between loose and stub flanges in a lap joint, see Figure 11.5-9 (typical); 
% B is inside diameter of flange; 
% B2 is the inside diameter of the contact face between loose and stub flanges in a lap joint, see Figure 11.5-9 (typical); 
% b is the effective gasket or joint seating width; 
% b0 is the basic gasket or joint seating width; 
% C is the bolt pitch circle diameter; 
% CF is the bolt pitch correction factor; 
% D is the inside diameter of shell; 
% db is bolt outside diameter; 
% dn is the bolt nominal diameter; 
% e is the minimum flange thickness, measured at the thinnest section; 
% fB is the bolt nominal design stress at operating temperature (see 11.4.3); 
% fB,A is the bolt nominal design stress at assembly temperature (see 11.4.3); 
% fH is the nominal design stress of the hub – see 11.5.4.2; 
% G is the diameter of gasket load reaction, as given by requirements in 11.5.2; 
% G1 is the assumed diameter of load reaction between loose and stub flanges in a lap joint; 
% g0 is the thickness of hub at small end; 
% g1 is the thickness of hub at back of flange; 
% H is the total hydrostatic end force; 
% HD is the hydrostatic end force applied via shell to flange; 
% HG is the compression load on gasket to ensure tight joint; 
% HT is the hydrostatic end force due to pressure on flange face; 
% h is the hub length; 
% hD is the radial distance from bolt circle to circle on which HD acts; 
% hG is the radial distance from gasket load reaction to bolt circle; 
% hL is the radial distance from bolt circle to circle on which load reaction acts for the loose flange in a 
% lap joint; 
% hT is the radial distance from bolt circle to circle on which HT acts; 
% K is the ratio of the flange diameters – see formulae 11.5-21 and 11.9-13; 
% k is stress factor defined in 11.5.4.2; 
% l0 is a length parameter given by Formula (11.5-22); 
% M is the moment exerted on the flange per unit of length, defined in 11.5.4.1; 
% MA is the total moment acting upon flange for assembly condition; 
% Mop is the total moment acting upon flange for operating condition; 
% m is a gasket factor; 
% Pe is the external calculation pressure, expressed as a positive number; 
% W is the design bolt load for assembly condition; 
% WA is the minimum required bolt load for assembly condition; 
% Wop is the minimum required bolt load for operating condition; 
% w is the contact width of gasket, as limited by gasket width and flange facing; 
% y is the minimum gasket or joint seating pressure; 
% betaF is a factor for integral method flange design as given in Figure 11.5-4; 
% betaFL is a factor for loose hubbed flanges as given in Figure 11.5-7; 
% betaT is a factor, given by formula (11.5-23); 
% betaU is a factor, given by formula (11.5-24); 
% betaV is a factor for the integral method, from Figure 11.5-5; 
% betaVL is a factor for loose hubbed flanges, from Figure 11.5-8; 
% betaY is a factor, given by Formula (11.5-25); 
% delta is the nominal gap between the shell and loose flange in a lap joint; 
% deltab is distance between centre lines of adjacent bolts; 
% lambda is a factor defined in 11.5.4.1; 
% sigmab is calculated bearing stress in a lap joint; 
% sigmaH is the calculated longitudinal stress in hub; 
% sigmar is the calculated radial stress in flange; 
% sigmatheta is the calculated tangential stress in flange; 
% phi is the hub stress correction factor for integral method flange design as given in Figure 11.5-6. 

MaterialS235JR;
fH = mat.fd;

bolt_table_ISO;
bolt.Rp02 = 580;
bolt.Rm = 830;

flange_type = 'stepped-bore'; % Fig 11.5-2
sealing_type = 'metal-to-metal';
% sealing_type = 'narrow-face-gasket';
% assembly_type = 'bolted-dome';
assembly_type = 'no-dome';


% Rp_bolts = min([Rp02_bolts/3, Rm_bolts/4]); Sec 11.4.3.1
% Rp_nuts = Rp02_nuts; Sec 11.4.3.1
fB = min([bolt.Rp02/3, bolt.Rm/4]);
fBA = fB;

f = fH; % see 11.5.4.2
nu = 0.3; %poisson



% B = 800; % inside diameter of flange
% C = 1010; % bolt pitch circle diameter
% A = 1060; % outside diameter of the flange
% gasket_diameter_mean = 880; % MY CHOICE

% %% TENTATIVE 1
% B = 1204; % inside diameter of flange
% C = 1300; % bolt pitch circle diameter
% A = 1375; % outside diameter of the flange
% gasket_diameter_mean = 1240; % MY CHOICE
% bolt_id_override = 30;


% %% SINGLE-CUP INTERNAL POSPRESSURE
% B = 1040; % inside diameter of flange
% gasket_diameter_mean = 1065; % MY CHOICE
% C = 1107.5; % bolt pitch circle diameter
% A = 1150; % outside diameter of the flange
% bolt_id_override = 23;
% n = 24; % number of bolts
% e = 20; % tentative thickness
% internal_pressure = true;
% P = 0.3; % pressure


%% SINGLE-CUP INTERNAL POSPRESSURE SMALL
B = 604; % inside diameter of flange
gasket_diameter_mean = 627.5; % MY CHOICE
C = 675; % bolt pitch circle diameter
A = 704; % outside diameter of the flange
bolt_id_override = 22;
n = 20; % number of bolts
e = 20; % tentative thickness
internal_pressure = true;
P = 0.3; % pressure

% %% SINGLE-CUP VACUUM
% B = 1208; % inside diameter of flange
% gasket_diameter_mean = 1240; % MY CHOICE
% C = 1300; % bolt pitch circle diameter
% A = 1350; % outside diameter of the flange
% bolt_id_override = 24;
% n = 24; % number of bolts
% e = 20; % tentative thickness
% internal_pressure = false;
% P = 0.1; % pressure



% % PN2.5
% B = 1234; % inside diameter of flange
% C = 1320; % bolt pitch circle diameter
% A = 1375; % outside diameter of the flange

% % PN6
% B = 1234; % inside diameter of flange
% C = 1340; % bolt pitch circle diameter
% A = 1405; % outside diameter of the flange

gasket_diameter_outer = gasket_diameter_mean + 4; % MY CHOICE
w = 20; % contact width of gasket (not used in metal-to-metal)
m = 0; % gasket factor
g0 = 4; % shell thickness
if matches(assembly_type, 'bolted-dome')
    a = 20; %  distance from top of flange to the mid-thickness line of the dome where it meets the flange
    R = 0.9*B;
end


if ~internal_pressure
    Mop_positive = 8.1031e+06; % see 11.5.5 Narrow face flanges subject to external pressure at final sentences
    % Mop_positive = 0;
end


if ~matches(sealing_type, 'metal-to-metal')
    y = 1.1; % the minimum gasket or joint seating pressure; see Annex H - Gasket factors m and y 
end


if matches(flange_type, 'stepped-bore') % WARNING: personal consideration
    g1 = 2*g0;
    h = g0; % weld height
end

%% Display

disp(['### Evaluating thickness: ', num2str(e), ' mm'])
if internal_pressure
    disp(['Internal pressure: ', num2str(P*10), ' bar'])
else
    disp(['External pressure: ', num2str(P*10), ' bar'])
end
disp(['Flange type: ', flange_type])


%% Bolt load and areas
b0 = w/2; %  eq (11.5-1) 
assert(b0>0.9, 'Units must be millimeters.')
if b0<=6.3 || matches(sealing_type, 'metal-to-metal')
    b = b0; 
    G = gasket_diameter_mean;
else
    b = 2.52*sqrt(b0);
    G = gasket_diameter_outer - 2*b;
end


%% Geometry
% hD = (C-B-g1)/2; 
assert(C>B)
assert(C>G)
hD = (C - B) / 2; % slip-on hubbed and stepped bore flanges (11.5-13)
hG = (C - G) / 2;
hT = (2*C-B-G)/4;

%% Loads
H = pi/4*(G^2*P);
HD = pi/4*(B^2*P);
HT = H - HD;

if matches(sealing_type, 'metal-to-metal')
    WA = 0;
    hR = (A - C)/2;
    MR = HD*hD + HT*hT;
    HR = MR/hR;
    Wop = H + HR;
    HG = 0;
else
    HR = 0;
    HG = 2*pi*G*b*m*P;
    WA = pi*b*G*y;
    Wop = H + HG;
end



ABmin = max([WA/fBA, Wop/fB]);
bolt_id = find(bolt_table.AB>(ABmin/n), 1, 'first');
if bolt_id>bolt_id_override
    error('Bolt override does not respect constraint.')
end
bolt_id = bolt_id_override; 
assert(~isempty(bolt_id), 'Cannot find bolt with proper dimensions');
assert(gasket_diameter_mean<C-bolt_table.dB0(bolt_id)/2, 'Cannot place the gasket')
AB = n*bolt_table.AB(bolt_id);
disp(['Bolt M', num2str(bolt_table.dB0(bolt_id)), ' [ID: ', num2str(bolt_id), ']'])
disp(['Bolt num: ', num2str(n)])

%% Flange Moments
if matches(sealing_type, 'metal-to-metal')
    e_min = sqrt(6*MR/(f*(pi*C-n*bolt_table.dh(bolt_id))));
    disp(['Minimum thickness by metal-to-metal formula: ', num2str(e_min), ' mm'])
    if e<e_min
        error(['Tentative thickness is too low. Minimum: ', num2str(e_min)])
    end
else
end

W = 0.5*(ABmin + AB)*fBA; % (11.5-16)

if matches(assembly_type, 'bolted-dome')
    Hr = HD*sqrt(4*R^2-B^2)/B;
    hr = e/2-a;
else
    Hr = 0;
    hr = 0;
end

% Assembly
% if matches(sealing_type, 'metal-to-metal')
%     MA = 0;
%     warning('Disabling assembly moment')
% else
    MA = W*hG; % assembly conditions  (11.5-17)
% end
if internal_pressure
    Mop = HD*hD + HT*hT + HG*hG - Hr*hr; % (11.5-18)
else
    Mop = max([HD*(hD - hG) + HT*(hT - hG), Mop_positive]);
end
disp(['MA (total moment, assembly): ', num2str(MA/1000), ' Nm'])
disp(['Mop (total moment, operating): ', num2str(Mop/1000), ' Nm'])

%% Flange stresses 11.5.4.1.1
deltab = pi*C/n;
db = bolt_table.dB0(bolt_id); % TO VERIFY
CF = max([1, sqrt(deltab/(2*db+6*e/(m + 0.5)))]); % (11.5-20)
K = A/B;

l0 = sqrt(B*g0);
betaT = (K^2*(1+8.55246*log10(K))-1)/((1.0472 + 1.9448*K^2)*(K-1));
betaU = (K^2*(1+8.55246*log10(K))-1)/(1.36136*(K^2-1)*(K-1));
betaY = (0.66845 + 5.7169*(K^2*log10(K))/(K^2-1))/(K-1);

MA_final = MA*CF/B;
Mop_final = Mop*CF/B;

load('flange_parameters/betaF.mat');
load('flange_parameters/betaV.mat');
load('flange_parameters/phi.mat');

betaF = betaF_struct.eval(g1/g0, h/l0, false);
betaV = betaV_struct.eval(g1/g0, h/l0, false);
phi = phi_struct.eval(g1/g0, h/l0, false);

disp(['betaF: ', num2str(betaF), ' | betaV: ', num2str(betaV), ' | phi: ', num2str(phi)]);


% [betaF, betaV, phi] = getFlangeParameters(g1, g0, h, l0, nu, false);
lambda = ((e*betaF + l0)/(betaT*l0) + (e^3*betaV)/(betaU*l0*g0^2));

sigmaH_fun = @(M) phi*M/lambda/g1^2;
sigmar_fun = @(M) (1.333*e*betaF+l0)*M/(lambda*e^2*l0);
sigmatheta_fun = @(M) betaY*M/e^2 - sigmar_fun(M)*(K^2+1)/(K^2-1);

sigmaH_A = sigmaH_fun(MA_final);
sigmar_A = sigmar_fun(MA_final);
sigmatheta_A = sigmatheta_fun(MA_final);

sigmaH_op = sigmaH_fun(Mop_final);
sigmar_op = sigmar_fun(Mop_final);
sigmatheta_op = sigmatheta_fun(Mop_final);

sigmaH_cases = [sigmaH_A, sigmaH_op];
sigmar_cases = [sigmar_A, sigmar_op];
sigmatheta_cases = [sigmatheta_A, sigmatheta_op];


name_cases = {'Assembly', 'Operating'};
for sel = 1:2

    sigmaH = sigmaH_cases(sel);
    sigmar = sigmar_cases(sel);
    sigmatheta = sigmatheta_cases(sel);
    name_case = name_cases(sel);

    if B<1000
        k = 1.0;
    elseif B<2000
        k = 2/3*(1+B/2000);
    else
        k = 1.333;
    end
    
    if sigmaH>1.5*min([f, fH])/k
        warning("NOT PASSED: Longitudinal stress, during " + name_case);
    end
    
    if sigmar>f/k
        warning("NOT PASSED: Radial stress, during " + name_case);
    end
    
    if sigmatheta>f/k
        warning("NOT PASSED: Tangential stress, during " + name_case);
    end
    
    if 0.5*(sigmaH + sigmar)>f/k
        warning("NOT PASSED: Longitudinal+Radial stress, during " + name_case);
    end
    
    if 0.5*(sigmaH + sigmatheta)>f/k
        warning("NOT PASSED: Longitudinal+Tangential stress, during " + name_case);
    end

end


%% Figures
if internal_pressure
    figure('NumberTitle', 'off', 'Name', 'Single Bolt Loads (W)'); hold on; grid on;
    W_labels = categorical({'Assembly','Operating'});
    W_values = [WA/n, 0; H/n, HG/n];
    W_bar = bar(W_labels, W_values, 'stacked');
    xtips1 = W_bar(end).XEndPoints;
    ytips1 = W_bar(end).YEndPoints;
    labels1 = string(W_bar(1).YData) + 'N';
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    ylabel('Bolt Loads [N]')
    xlabel('Load Condition')
end


figure('NumberTitle', 'off', 'Name', 'Parameters');
subplot(3,1,1);
betaF_ax = betaF_struct.plot();
plot(betaF_ax, g1/g0, betaF, 'ko');
title('BetaF')
subplot(3,1,2);
betaV_ax = betaV_struct.plot();
plot(betaV_ax, g1/g0, betaV, 'ko');
title('BetaV')
subplot(3,1,3);
phi_ax = phi_struct.plot();
plot(phi_ax, g1/g0, phi, 'ko');
title('phi')


figure('NumberTitle', 'off', 'Name', 'Stresses'); hold on; grid on;
stress_limits = [1.5*min([f, fH])/k, f/k, f/k, f/k, f/k];
plot(stress_limits, '-ro');
plot([sigmaH_A, sigmar_A, sigmatheta_A, 0.5*(sigmaH_A + sigmar_A), 0.5*(sigmaH_A + sigmatheta_A)], '-bx');
plot([sigmaH_op, sigmar_op, sigmatheta_op, 0.5*(sigmaH_op + sigmar_op), 0.5*(sigmaH_op + sigmatheta_op)], '-gx');
ylabel('Stress [MPa]')
xlabel('Load Condition')
set(gca, 'XTick',1:numel(stress_limits), 'xticklabel',{'Long'; 'Radial'; 'Tang'; '0.5*(Long+Radial)'; '0.5*(Long+Tang)'})
legend('LIMITS', 'Assembly', 'Operational')

force_scale = 1/500/n;
eG = 3; % gasket thickness
figure('NumberTitle', 'off', 'Name', 'Sketch'); axis equal; hold on;
rectangle('Position', [-A/2, 0, A/2-B/2, e])
rectangle('Position', [-C/2 - db/2, 0, db, e])
if matches(sealing_type, 'metal-to-metal')
    rectangle('Position', [-G/2-eG/2, -eG/2, eG, eG], 'Curvature', [1 1], 'FaceColor','k')
    line([-A/2, -B/2], [0, 0], 'LineWidth', 2)
else
    rectangle('Position', [-G/2 - eG/2, -eG, w, eG], 'FaceColor','b')
end
if matches(flange_type, 'stepped-bore')
    rectangle('Position', [-B/2, e/2, g0, e])
    patch([-B/2-h, -B/2, -B/2], [e, e+h, e], [0 0 0])
    patch([-B/2, -B/2+h, -B/2], [e/2, e/2, e/2-h], [0 0 0])
end
line([-C/2, -C/2], [0, e],'LineStyle','-.', 'Color', 'k')
quiver(-C/2, e+W*force_scale, 0, -W, force_scale);
quiver(-C/2+hG, -HG*force_scale, 0, HG, force_scale)
quiver(-C/2+hT, -HT*force_scale, 0, HT, force_scale)
quiver(-B/2*0.95, -H*force_scale, 0, H, force_scale)







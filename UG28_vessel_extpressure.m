clear
close all
clc

unit_conversion;

%% Parameters
% A  = 1.0; % factor “A,” strain, from ASME Section II, Part D, Subpart 3, dimensionless
% As  = 1.0; % cross-sectional area of stiffener, in.2
% B  = 1.0; % factor “B,” allowable compressive stress, from ASME Section II, Part D, Subpart 3, psi
% DL  = 1.0; % outside diameter of the large end of cone, in.
% Ds  = 1.0; % outside diameter of small end of cone, in.
% 
% I  = 1.0; % actual moment of inertia of stiffener, in.4
% Is  = 1.0; % required moment of inertia of stiffener, in.4
% I0s  = 1.0; % required moment of inertia of combined shell-ring cross section, in.4
% L  = 1.0; % for cylinders  - the design length for external pressure, including 1/3 the depth of heads, in. For cones - the design length for external pressure (see Figures 2-1b and 2-1c), in.
% Le  = 1.0; % equivalent length of conical section, in.
% Ls  = 1.0; % length between stiffeners, in.
% P  = 1.0; % design internal pressure, psi
% Pa  = 1.0; % allowable external pressure, psi
% Px  = 1.0; % design external pressure, psi
% Ro  = 1.0; % outside radius of spheres and hemispheres, crown radius of torispherical heads, in.
% % t  = 1.0; 
% te  = 1.0; % equivalent thickness of cone, in.
% a  = 1.0; % half apex angle of cone, degrees

%% Set parameters
E  = 210e9 * FROM_PA_TO_PSI; % modulus of elasticity, psi
D  = 1200e-3 * FROM_M_TO_INCHES; % inside diameter of cylinder, in.
LTT  = 0.8 * FROM_M_TO_INCHES; % length of straight portion of shell, tangent to tangent, in.
t = 4e-3 * FROM_M_TO_INCHES; % thickness of cylinder, head or conical section, in.
h = 192e-3 * FROM_M_TO_INCHES; % total the depth of the head from the tangent line (see Figure UG-28)

%% Calculated parameters
Do  = D + 2*t; % outside diameter of cylinder, in.

if Do/t<10
    error('Unsupported case')
end

%% Step 1-3
% L = LTT + D/3; % hemi-heads
% L = LTT + D/6; % 2:1 S.E. heads
% L = LTT + 0.112*D; % 100% e 6% heads
L = LTT + h/3;
% 
% disp(['L/Do: ', num2str(L/Do)])
% disp(['Do/t: ', num2str(Do/t)])

%% Step 4
% from Fig. 2-1e  Geometric chart for components under external or compressive loadings (for all materials). (ASME Code, Section VIII, Div. 1.)
A = getFactorA(Do/t, L/Do, true); % factor "A," strain, from ASME Section II, Part D, Subpart 3, dimensionless
[B, applicableB] = getFactorB(A, true);

%% Step 5
if applicableB
    Pa = 4*B/3/(Do/t); % allowable external pressure, psi
else
    Pa = 2*A*E/3/(Do/t); % allowable external pressure, psi
end

disp(['Allowable external pressure: ', num2str(Pa/FROM_PA_TO_PSI*1e-5), ' bar'])


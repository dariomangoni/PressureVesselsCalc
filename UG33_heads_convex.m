clear
close all
clc

unit_conversion;

% Option 1 - use UG32 with a pressure 1.67 higher
% Option 2 - use formulas for ellipsoidal heads with appropriate values of Ro

% Ro = for hemispherical heads, the outside radius
% = for ellipsoidal heads, the equivalent outside spherical radius taken as KoDo
% = for torispherical heads, the outside radius of the crown portion of the head
% Do = outside diameter of the head skirt
% Pa: maximum allowable external working pressure
t = 4e-3 * FROM_M_TO_INCHES;
Do = 1200e-3 * FROM_M_TO_INCHES + 2*t;
Ro = 0.9 * Do;
L = Ro;

% A = 0.125/(Ro/t); why this?
disp(['L/Do: ', num2str(L/Do)])
disp(['Do/t: ', num2str(Do/t)])
A = getFactorA(Do/t, L/Do, true);
[B, applicableB] = getFactorB(A, true);

if applicableB
    Pa = B/(Ro/t);
else
    Pa = 0.0625*E/(Ro/t).^2;
end
disp(['Pa: ', num2str(Pa/FROM_PA_TO_PSI*1e-5), ' bar'])

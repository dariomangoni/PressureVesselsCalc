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

assert(e/(Di+2*e)<=0.16, 'Sec 7.4.1 - Formula not valid for this sizes.')

ea = e;
Dm = Di + e;

P = 2*f*z*ea/Dm;

disp(['Allowable pressure: ', num2str(P*10), ' bar'])

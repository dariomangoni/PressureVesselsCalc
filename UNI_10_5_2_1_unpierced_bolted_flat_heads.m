clear
close all
clc

bolt_table_ISO;
bolt_id = 28; % (taken from UNI_flange.m)

% tB is the mean bolt pitch in a bolted flat end; WARNING: it seems that this is not the thread pitch, but is the distance between bolts (see variable C)

MaterialS235JR; % TO CHANGE WITH PROPER MATERIAL

f = mat.fd;
fA = mat.fd;
G = 1220; % gasket reaction diameter
B = 1200;
C = 1300; % bolt pitch circle
A = 1400;
b = 6;
P = 6/10;
m = 0;
nu = 0.3;

W = 1; % TODO
fA = 1; % TODO

e1a = 22; % analysis thickness for the flanged extension

n = 20; % number of bolts

tB = C/n; % Bs in ASME
db = bolt_table.dB0(bolt_id); % assuming db = dB0; since outer is nominal, isn't it?

%% thickness within gasket shell
CF = max([1, sqrt(tB/(2*db+6*e1a/(m + 0.5)))]);
eA = sqrt(CF*3*(C-G)/(pi*G)*W/fA);
eP = sqrt((3*(3+nu)/32*G^2 + 3*CF*(G/4 + 2*b*m)*(C-G))*P/f);
e = max([eA, eP]);


%% thickness for the flanged extension
eP1 = sqrt(3*CF*(G/4+2*b*m)*(C-G)*P/f);
e1 = max([eA, eP1]);

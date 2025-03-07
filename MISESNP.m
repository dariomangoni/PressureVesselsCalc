
clear
close all
clc

L = 0.8; % unsupported length of cylinder
D = 1.2;  % mean shell diameter
A = D/2; % mean shell radius
H = 3e-3; % wall thickness of shell
E = 210e9; % Young2s modulus of elasticity
NU = 0.3; %Poisson0s ratio
SIGMAYP = 300e6; % yield stress

C1 = (L / D) ^ 2;
C2 = (H / D) ^ 3;
C3 = SIGMAYP / E;
LAMBDA = SQR(SQR(C1 / C2)) * SQR(C3);
C1 = (H / D) ^ 2.5;
C2 = SQR(H / D);
C3 = L / D;
P = 2.6 * E * C1 / (C3 - .45 * C2);
disp(['BUCKLING PRESSURE(DTMB): ', num2str(P)]);
P = 1E+20;
for N = 1:200
    N2 = N * N;
    C1 = pi * A / L;
    NU2 = (1 - NU ^ 2);
    C2 = (N2 - 1 + C1 ^ 2) ^ 2;
    C2 = H ^ 2 * C2 / (12 * A ^ 2 * NU2);
    C3 = (N2 * (1 / C1) ^ 2 + 1) ^ 2;
    C4 = N2 - 1 + .5 * C1 ^ 2;
    PCR = E * H / (A * C4) * (1 / C3 + C2);
    if P < PCR
        disp(['BUCKLING PRESSURE(VON MISES): ', num2str(P), ' LOBES: ', num2str(N - 1)])
        break
    end
    P = PCR;
end
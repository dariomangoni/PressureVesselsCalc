clear
close all
clc

bolt_table_ISO;

flange_ring_type = 'integral';
% flange_type = 'blank';
% flange_type = 'loose';

% connected_shell = 'tapered-hub'; % 5.2.3.1
connected_shell = 'no-hub'; % 5.2.3.2
% connected_shell = 'blank'; % 5.2.3.3
% connected_shell = 'collar-stub'; % 5.2.3.4


nB = 36;

%% UNI 1591-3 - 4.2.2 Flange ring
% 4.2.2.1 Bolt holes 
pB = pi*d3/nB; % pitch between bolts
d5e = d5 * sqrt(d5/pB); % Effective diameter of bolt hole: 
d5 = d5t * l5t/eFb; % Diameter of blind holes is assumed to be: 
d3e = d3 * (1-2/nB^2); % Effective bolt circle diameter:

% 4.2.2.2 Effective dimensions of flange ring 
if matches(flange_ring_type, 'integral') || matches(flange_ring_type, 'blank')
    bF = (d4-d0)/2 - d5e;
    bL = 0;
    dL = 0;
    eL = 0;
    dF = (d4+d0)/2;
    eF = 2*AF/(d4-d0);
elseif matches(flange_ring_type, 'loose')
end

% 4.2.3 Connected shell 
if matches(connected_shell, 'no-hub')
    % 4.2.3.2 Flange without hub 
    eE = eS;
    dE = dS;
elseif matches(connected_shell, 'blank')
    eE = 0;
    dE = d0;
else
    error('Case not supported')
end

%% UNI 1591-3 - 4.2.4 Lever arms 
hP = ((dGe - dE)^2*(2*dGe+dE)/6+2*ep^2*dF)/dGe^2; % (13, pag 19 PDF)
hD = (dMe - dGe)/2;

if matches(flange_ring_type, 'integral') || matches(flange_ring_type, 'blank')
    % 4.2.4.2 Integral flange and blank flange
    hG = (d3e - dGe)/2;
    hH = (d3e + dE)/2;
    hL = 0;
    hM = (d3e - dMe)/2;
elseif matches(flange_ring_type, 'loose')
    error('Case not supported')
end

%% 4.2.5 Flexibility-related flange parameters
% NOT APPLICABLE

%% 4.3 Bolt parameters 

% UNI 1591-1 - 5.3.2 | UNI 1591-3 - 4.3.2 Effective cross-section area of bolts
bolt_id = 4;
dB0 = bolt_table.dB0(bolt_id);
dBe = bolt_table.dBe(bolt_id);
dBs = bolt_table.dBs(bolt_id);
AB = bolt_table.AB(bolt_id);

XB = 4*(ls/dBs^2 + le/dBe^2 + 0.8/dB0)/(nB*pi);

%% UNI1591-3 4.4 Gasket parameters | UNI13445-3 Annex G.5.3
bGt = (dG2 - dG1)/2;
dGt = (dG2 + dG1)/2;
AGt = pi*dGt*bGt;





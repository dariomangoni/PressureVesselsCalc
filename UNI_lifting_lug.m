clear
close all
clc


a1 = 15; % eccentricity of load; 
a2 = 20; % distance from load to shell or reinforcing plate; 
b1 = 80; % length of lifting lug at shell junction; 
b2 = 80; % width of reinforcing plate; 
b3 = 100; % length of reinforcing plate; 
% x = ; % distance between the axis of semi-ellipsoidal head and the centre of the lifting lug ; 
% FR = ; % local force on a shell; 
% FRmax = ; % maximum allowable local force on a shell; 
Gmax = 15000; % total vessel weight [N]; 
% beta = 135*pi/180; % angle between direction of force and normal to the shell; 
e2 = 6; % reinforcement thickness

lifting_joint_point_height = 4000;
lifting_lugs_num = 3;

ea = 4; % shell thickness
en = 4; % shell nominal thickness
Di = 1200; % internal diameter of shell
Deq = Di; % for cylindrical shells
P = 0; % pressure (pos if internal, negative if external)
MaterialS235JR;

f_test = mat.ftest;
f = mat.fd;

F = Gmax;

direction = 'longitudinal';
reinforcing_plate = true;
% condition = 'design';
condition = 'test|transport|lifting';

%% APPLICATION-SPECIFIC DATA
if ~reinforcing_plate
    e2 = 0;
end
lifting_lugs_hole_distance = Di + 2*ea + 2*e2 + 2*a2;

beta = atan2(lifting_joint_point_height, -lifting_lugs_hole_distance/2);

FR = Gmax/sin(beta)/lifting_lugs_num;

F = FR*sin(beta);
M = F*(a2+e2);

%% CALCULATIONS
if reinforcing_plate
    b = b3;
    assert(e2>=en)
    assert(b3<=1.5*b1)
    assert(b3>=b1+2*e2) % my consideration
    K15 = min([]);
else
    b = b1;
    assert(b3>=b1+2*ea) % my consideration
end


assert(0.001<=en/Deq)
assert(en/Deq<=0.05)

sigma_mx_plus = P*Deq/4/ea + 1/(pi*Deq*ea)*(F+4*M/Deq); % (16.6-9)
sigma_mx_mins = P*Deq/4/ea + 1/(pi*Deq*ea)*(F-4*M/Deq); % (16.6-9)
sigma_my = P*Deq/2/ea;

if matches(condition, 'design')
    K2 = 1.25;
elseif matches(condition, 'test|transport|lifting')
    K2 = 1.05;
    f = f_test;
else
    error('condition not covered')
end

if matches(direction, 'longitudinal')
    lambda1 = b/sqrt(Deq*ea);
    lambda = lambda1;
    nu1 = min([0.08*lambda1, 0.2]);
    sigma_m = sigma_my;
    nu2 = sigma_m/K2/f; % eq (16.6-8) 
    K1 = (1-nu2^2)/(1/3+nu1*nu2+sqrt((1/3+nu1*nu2)^2 + (1-nu2^2)*nu1^2)); % eq (16.6-7) 
    K13 = 1/1.2/sqrt(1+0.06*lambda^2);
    K14 = 1/0.6/sqrt(1+0.03*lambda^2);
    K15 = min([2.0, 1+2.60*(Deq/ea)^0.30*b2/Deq]); % (16.7-2) 
elseif matches(direction, 'circumferential')
    K15 = min([1.8, 1+2.65*(Deq/ea)^0.33*b2/Deq]); % (16.7-3)
    error('Case not included')
end

sigmab_all = K1*K2*f; 

FLmax = sigmab_all * ea^2/K13; % (16.6-21)
MLmax = sigmab_all * ea^2*b/K14; % (16.6-22)

if reinforcing_plate
    FRmax = K15*sigmab_all*ea^2/(K13*abs(cos(beta))+K14*abs((a2+e2)*sin(beta)-a1*cos(beta))/b3);
else
    FRmax = sigmab_all*ea^2/(K13*abs(cos(beta))+K14*abs(a2*sin(beta)-a1*cos(beta))/b1);
end

figure('NumberTitle', 'off', 'Name', 'FR');
bar(["FR", "FRmax"],[FR, FRmax])

disp(['FR: ', num2str(FR), ' N (max: ', num2str(FRmax), ' N)'])
assert(FR<FRmax)


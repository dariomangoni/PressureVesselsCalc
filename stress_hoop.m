clear
close all
clc

%% general parameters
t = 3e-3; % thickness
p = 3e5; % pressure (relative), positive if internal

%% thin-walled parameters
r = 1200e-3; % radius (inner radius if p positive, outer radius negative otherwise)

%% thick-walled parameters

r1 = 1200e-3; % inner radius
r2 = 1203e-3; % outer radius

if t/r <= 1/30
    % Thin-walled circular cylinder
    stress_hoop_thin_cylinder = p*r/t; % [Ross, pag. 15]
    stress_long_thin_cylinder = p*r/2/t; % [Ross, pag. 15]
    
    % Thin-walled sphere
    stress_hoop_thin_sphere = p*r/2/t; % [Ross, pag. 16]
    stress_long_thin_sphere = p*r/2/t; % [Ross, pag. 16]
else

    % Thick-walled circular cylinder (max stress at the internal radius)
    if p>0
        stress_hoop_thick_cylinder = p*(r1^2 + r2^2)/(r2^2 - r1^2); % [Ross, pag. 16]
    else
        stress_hoop_thick_cylinder = -2*p*(r2^2)/(r2^2 - r1^2); % [Ross, pag. 16]

    end

end

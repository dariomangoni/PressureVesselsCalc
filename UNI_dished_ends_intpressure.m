clear
close all
clc

% De: is the outside diameter of the cylindrical flange; 
% Di: is the inside diameter of the cylindrical flange; 
% eb: is required thickness of knuckle to avoid plastic buckling; 
% es: is required thickness of end to limit membrane stress in central part; 
% ey: is required thickness of knuckle to avoid axisymmetric yielding; 
% fb: is design stress for buckling formula; 
% hi: is inside height of end measured from the tangent line; 
% K: is shape factor for an ellipsoidal end as defined in Formula (7.5-18); 
% N: is a parameter defined by Formula (7.5-12); 
% R: is inside spherical radius of central part of torispherical end; 
% X: is ratio of knuckle inside radius to shell inside diameter; 
% Y: is a parameter defined by Formula (7.5-9); 
% Z: is a parameter defined by Formula (7.5-10); 
% beta: is a factor given by Figures 7.5-1 and 7.5-2 or by the procedure in 7.5.3.5.
% r: is the inside radius of curvature of a knuckle. 
% Rp02T = 355; % MPa

% f: nominal design stress

%% Parameters
MaterialS235JR;

mode = 'oper';
Di = 1200;
R = Di;
r = 75; % TO VERIFY
P = 1/10; % pressure
f = fd;
z = 0.85; % joint coefficient Table 5.6-1

internal_pressure = true;


assert(r <= 0.2*Di)
assert(r >= 0.06*Di)


for mode_sel = 1:2
    e = 4; % reset

    for iter = 1:20
    
    
        De = Di + 2*e;
        
    
        assert(r >= 2*e)
        assert(e <= 0.08*De)
        assert(e >= 0.001*De)
        assert(R >= 2*e)
        
        %% Design
        
        X = r/Di;
        Y = min([e/R, 0.04]);
        

        if internal_pressure
            N = 1.006 - 1/(6.2 + (90*Y)^4);
        else
            N = 1; % 8.8.2 Torispherical ends 
        end

        Z = log10(1/Y);
        beta006 = N*(-0.3635*Z^3+2.2124*Z^2-3.2937*Z+1.8873);
        beta010 = N*(-0.1833*Z^3+1.0383*Z^2-1.2943*Z+0.837);
        beta020 = max([0.95*(0.56-1.94*Y-82.5*Y^2), 0.5]);
        if X < 0.06
            assert(false, 'X cannot be smaller than 0.06')
        elseif X <= 0.06+1e-3
            beta = beta006;
        elseif X < 0.1
            beta = 25*((0.1-X)*beta006 + (X-0.06)*beta010);
        elseif X < 0.1+1e-3
            beta = beta010;
        elseif X < 0.2
            beta = 10*((0.2-X)*beta010 + (X-0.1)*beta020);
        elseif X <0.2+1e-3
            beta = beta020;
        else
            assert(false, 'X cannot be greater than 0.2')
        end

        
        es = P*R/(2*f*z - 0.5*P);
        ey = beta*P*(0.75*R + 0.2*Di)/f;
        if mode_sel == 1 % operating conditions
            fb = Rp02T/1.5;
        elseif mode_sel == 2 % testing conditions
            fb = Rp02T/1.05;
        else
            error('Plese set ''mode'' appropriately')
        end
        eb = (0.75*R+0.2*Di) * (P/(111*fb)*(Di/r)^0.825)^(1/1.5);
        
        e_new = max([es, ey, eb]);
    
    
        if abs(e_new-e)<0.001
            if mode_sel == 1 % operating conditions
                disp(['Required thickness for pressure in operating conditions of ', num2str(P*10), ' bar: ', num2str(e), ' mm'])
            elseif mode_sel == 2 % testing conditions
                disp(['Required thickness for pressure in testing conditions of ', num2str(P*10), ' bar: ', num2str(e), ' mm'])
            end
            break
        end
        e = e_new;
    
    end
end

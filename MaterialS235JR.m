

mat.Rp02T = 235; % MPa, 0,2% proof strength at calculation temperature
mat.Rm20 = 400; % tensile strength

mat.fd = min([mat.Rp02T/1.5, mat.Rm20/2.4]);
mat.ftest = mat.Rp02T/1.05;

if exist('f', 'var')
    if f > mat.ftest && f>mat.fd
        error('Selected f is above material thresholds.')
    else
        if f<mat.ftest
            disp('Available for testing conditions.')
        end
        if f<mat.fd
            disp('Available for operating conditions.')
        end
    end
end

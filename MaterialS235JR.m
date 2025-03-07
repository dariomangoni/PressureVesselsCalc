

Rp02T = 235; % MPa, 0,2% proof strength at calculation temperature
Rm20 = 400; % tensile strength

fd = min([Rp02T/1.5, Rm20/2.4]);
ftest = Rp02T/1.05;

if exist('f', 'var')
    if f > ftest && f>fd
        error('Selected f is above material thresholds.')
    else
        if f<ftest
            disp('Available for testing conditions.')
        end
        if f<fd
            disp('Available for operating conditions.')
        end
    end
end

clear
close all
clc


betaF_struct = getDataFromWebplotDigitizer('betaF.csv');
betaV_struct = getDataFromWebplotDigitizer('betaV.csv');
phi_struct = getDataFromWebplotDigitizer('phi.csv');

betaF_struct.headers = betaF_struct.headers/100;
betaV_struct.headers = betaV_struct.headers/100;
phi_struct.headers = phi_struct.headers/100;

betaF_struct.eval = @(x, par, fail_if_param_out) evalParametricFunction(betaF_struct, x, par, fail_if_param_out);
betaV_struct.eval = @(x, par, fail_if_param_out) evalParametricFunction(betaV_struct, x, par, fail_if_param_out);
phi_struct.eval = @(x, par, fail_if_param_out) evalParametricFunction(phi_struct, x, par, fail_if_param_out);


betaF_struct.plot = @() plotWebplotdigitizer(betaF_struct, @(x, y) semilogx(x, y, '-x'));
betaV_struct.plot = @() plotWebplotdigitizer(betaV_struct, @(x, y) semilogx(x, y, '-x'));
phi_struct.plot = @() plotWebplotdigitizer(phi_struct, @(x, y) loglog(x, y, '-x'));

figure;
betaF_struct.plot();
figure;
betaV_struct.plot();
figure;
phi_struct.plot();

save('betaF.mat', 'betaF_struct')
save('betaV.mat', 'betaV_struct')
save('phi.mat', 'phi_struct')


function ax = plotWebplotdigitizer(val_struct, plotfun)
    leg = cell(1, numel(val_struct.headers));
    for sel = 1:numel(val_struct.headers)
        plotfun(val_struct.values{sel}(:,1), val_struct.values{sel}(:,2));
        hold on; grid on;
        leg{sel} = num2str(val_struct.headers(sel));
    end
    legend(leg);
    ax = gca;
end




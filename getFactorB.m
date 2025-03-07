function [B, applicable] = getFactorB(A, plotFlag)


if nargin<2
    plotFlag = false;
end


tableB = ...
[0.135e-04, 0.200e+03;
0.645e-03, 0.940e+04;
0.700e-03, 0.990e+04;
0.800e-03, 0.105e+05;
0.900e-03, 0.109e+05;
0.100e-02, 0.112e+05;
0.150e-02, 0.119e+05;
0.200e-02, 0.123e+05;
0.950e-02, 0.138e+05;
0.200e-01, 0.138e+05;
0.100, 0.138e+05];


applicable = true;
if A<min(tableB(:,1))
    applicable = false;
    B = NaN;
    return;
elseif A>=max(tableB(:,1))
    B = tableB(end,2);
else
    B = interp1(tableB(:,1), tableB(:,2), A);
end

if plotFlag
    figure('NumberTitle', 'off', 'Name', 'Factor B');
    loglog(tableB(:,1), tableB(:,2), 'Marker', 'x');
    hold on; grid on;
    loglog(A, B, 'ko');
    legend({'<300 °F = 150 °C', 'Picked point'})
    ylabel('B')
    xlabel('A')
end

end
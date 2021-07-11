arg = linspace(-1,1);
arg_num = length(arg);
max_degree = 5;

PD_exact = zeros(max_degree+1, arg_num);
PD_exact(1,:) = 0;
PD_exact(2,:) = 1;
PD_exact(3,:) = 3 .* arg;
PD_exact(4,:) = 1/2 .* (15*arg.^2 - 3);
PD_exact(5,:) = 1/8 .* (4*35*arg.^3-60*arg);
PD_exact(6,:) = 1/8 .* (5*63*arg.^4-3*70*arg.^2+15);

PD_cal = LegendrePolynomialD(max_degree, arg);

relerr = abs(PD_exact - PD_cal) ./ abs(PD_exact + PD_cal);

%% compare the results
fig = Figure;
for i = 1:max_degree+1
    subplot(max_degree+1,1,i);
    plot(arg, PD_exact(i,:));
    hold on
    plot(arg, PD_cal(i,:), ':');
    fig.Init;
end

fig_relerr = Figure;
plot(arg, relerr.');
fig_relerr.Init;
title('Relative error');

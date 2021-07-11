arg = linspace(-1,1);
arg_num = length(arg);
max_degree = 5;

PD_exact = LegendrePolynomialDClosdeForm(max_degree, arg);

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

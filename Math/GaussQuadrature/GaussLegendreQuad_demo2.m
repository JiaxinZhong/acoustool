% Example: integral of exp(-x) from a to infinity
gauss_num = (3:2e1).';
a = 1;
res = 0 * gauss_num;
for i = 1:length(gauss_num)
    gauss = GaussLegendreQuadParam(gauss_num(i), 'lower_limit', a, ...
        'is_singular', 1);

    integrand = exp(-gauss.zero);
    res(i) = sum(integrand .* gauss.weight);
end
fig = Figure;
plot(gauss_num, res);

res_correct = exp(-a);
    
fprintf("Correct: %e\n", res_correct);
fprintf("Calculated: %e\n", res(end));

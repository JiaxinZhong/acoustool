% Example: integral of exp(-x^2/2) from -1 to 1
%   see Example 5.15 in page 276 of Numerical Analysis 2ed
%   written by Timothy Sauer

gauss_num = (3:20).';
res = 0 * gauss_num;
for i = 1:length(gauss_num)
    gauss = GaussLegendreQuadParam(gauss_num(i), ...
        'lower_limit', -1, 'upper_limit', 1);

    integrand = exp(-gauss.zero.^2 / 2);
    res(i) = sum(integrand .* gauss.weight);
    
end
fig = Figure;
plot(gauss_num, res);

res_correct = 1.71124878378430;


fprintf("Correct: %e\n", res_correct);
fprintf("Calculated: %e\n", res(end));

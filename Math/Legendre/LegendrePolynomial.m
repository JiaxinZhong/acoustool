% ==============================================================================
% Calculate Legendre polynomials
% ------------------------------------------------------------------------------
% INPUT
%	n			--- degree; must be in dimension 1
%	x			---	argument; must be in dimensions 2, 3, ...
% ------------------------------------------------------------------------------
% OUTPUT
%	pn			--- Legendre polynomial
% ==============================================================================
function pn = LegendrePolynomial(max_degree, x)

% 	if ~iscolumn(n)
% 		error('n必须在第1维！\n');
% 	end
	if (~isscalar(x) && size(x, 1) > 1)
		error('宗量必须在第2维及以上！\n')
	end

	% 先算出0到N_MAX度的值，再取目标值n
	x_row = x(:).';

	% n为第1维，x为第2维
% 	pn = zeros(N_MAX+1, length(x));
    degree = (0:max_degree).';
    pn = 0 * degree .* x_row;
	% pnDerivative = zeros(n+1,1);

	% 宗量接近于1时，特殊处理
	% if abs(x-1)<10*eps
		% pn = ones(N_MAX+1,1);
		% n_span = (0:N_MAX).';
		% % pnDerivative = 1/2*n_span.*(n_span+1);
		% % legendrePolynomial.data = pn;
		% % legendrePolynomial.derivative = pnDerivative;
		% return
	% end

	% 初始值
	pn(1, :) = 1;
    if (max_degree > 0)
        pn(2, :) = x_row;

	for nn = 2:max_degree
		pn(nn+1,:) = ((2*nn-1).*x_row.*pn(nn-1+1,:)-(nn-1)*pn(nn-2+1,:))/nn;
    end

    
%     pn = pn(max_degree+1, :);
	pn = reshape(pn, size(0*degree.*x));
end

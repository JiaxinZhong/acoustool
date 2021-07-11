% ==============================================================================
% Calculate the derivative of the Legendre polynomial 
% ------------------------------------------------------------------------------
% INPUT
%	max_degree  --- The maximum degree
%	arg         --- The argument. 
%                   Requirements: 1. size(arg, 1) == 1
% ------------------------------------------------------------------------------
% OUTPUT
%   PD          --- The result with the size of (max_degree+1) x dim(arg).
%                   The column corresponds to the degree 0, 1, 2, ..., max_degree
% ==============================================================================
function PD = LegendrePolynomialD(max_degree, arg)

	if size(arg, 1) > 1
		error('Requirement not satisfied: size(arg,1) == 1\n')
	end

	arg_row = arg(:).';
	leg_n = LegendrePolynomial(max_degree+1, arg_row);

	degree = (0:max_degree).';
    % valid only when abs(arg) != 1
	PD = (degree+1) ./ (1-arg_row.^2) .* ...
        (arg_row .* leg_n(degree+1,:) - leg_n(degree+2,:));
    
    % special case
    special_case = abs(arg_row) == 1;
    num = sum(special_case);
    PD(:, special_case) = 1/2 .* (arg_row(special_case)).^(degree+1) .* degree .* (degree + 1) + zeros(1,num);

	PD = reshape(PD, size(degree.*arg));
end

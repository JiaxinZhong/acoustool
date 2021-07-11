% ==============================================================================
% Calculate the derivative of the Legendre polynomial
% Using the closed form method
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
function PD = LegendrePolynomialDClosdeForm(max_degree, arg)

    MAX_DEGREE = 5;

    if (max_degree > MAX_DEGREE)
        error('The maximum degree should be no larger than %d.\n', MAX_DEGREE);
    end
	if size(arg, 1) > 1
		error('Requirement not satisfied: size(arg,1) == 1\n')
	end

    degree = (0:MAX_DEGREE).';
	arg_row = arg(:).';
    PD = 0 * repmat(arg_row, MAX_DEGREE+1, 1);

    PD(1,:) = 0;
    PD(2,:) = 1;
    PD(3,:) = 3 .* arg_row;
    PD(4,:) = 1/2 .* (15*arg_row.^2 - 3);
    PD(5,:) = 1/8 .* (4*35*arg_row.^3-60*arg_row);
    PD(6,:) = 1/8 .* (5*63*arg_row.^4-3*70*arg_row.^2+15);

	PD = reshape(PD, size(degree.*arg));
end

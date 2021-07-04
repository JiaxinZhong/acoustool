% =========================================================================
% Import the parameters (zeros and weights) for Gauss-Legendre Quadrature
%	method.
% -------------------------------------------------------------------------
% INPUT
%	GaussNumber	---	scalar
%   lower_limit, upper_limit
%               --- The range of the integral.
%   Dimension:  gauss => lower_limit .* upper_limit
% -------------------------------------------------------------------------
% Output
%	zero, weight	--- zeros and weights
% =========================================================================

function param = GaussLegendreQuadParam(gauss_num, varargin)

    p = inputParser();
    % two number or one number
    % The range is -1 to 1 if the range of the integral is not assigned
    addParameter(p, 'lower_limit', -1);
    addParameter(p, 'upper_limit', 1);
    addParameter(p, 'is_singular', 0);
    parse(p, varargin{:});
    ip = p.Results;

    if size(ip.lower_limit,1) > 1
        error('Must: size(ip.lower_limit,1) > 1)\n');
    end
    if size(ip.upper_limit,1) > 1
        error('Must: size(ip.upper_limit,1) > 1)\n');
    end
    
	if numel(gauss_num) > 1
		error('numel(GaussNumber) must be 1!\n')
    end

    if gauss_num > 1e3
        gauss_num = 1e3;
        fprintf('The input of Gauss number is too large! It is set as 1000!\n')
    end
    
	fn_zero = sprintf('GaussLegendreQuadrature_zero_%s.txt', ...
		num2str(gauss_num));
	fn_weight = sprintf('GaussLegendreQuadrature_weight_%s.txt',...
		num2str(gauss_num));
    param.zero = importdata(fn_zero);
	param.zero = param.zero(:);
    param.weight = importdata(fn_weight);
	param.weight = param.weight(:);

    switch ip.is_singular
        case 0
            param.weight = param.weight .* (ip.upper_limit - ip.lower_limit)/2;
            param.zero = ((ip.upper_limit - ip.lower_limit) ...
                .* param.zero + (ip.upper_limit + ip.lower_limit))/2;
        case 1
            % When the range is inifnite, ip.upper_limit is invalid.
            % Using trigometric transformations
            param.weight = param.weight .* pi/4 .* (sec(pi/4*(param.zero + 1))).^2;
            param.zero = tan(pi/4 * (param.zero + 1)) + ip.lower_limit;
    end
end

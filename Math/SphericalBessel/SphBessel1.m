% =========================================================================
% Calculate the spherical Bessel function of the first kind
% -------------------------------------------------------------------------
% INPUT
%	n	--- order
%	arg	--- argument
% -------------------------------------------------------------------------
% OUTPUT
%	j	--- the results
% =========================================================================
function j = SphBessel1(n, arg, varargin)
	
    p = inputParser;
    addParameter(p, 'Method', 'Normaliarged');
    addParameter(p, 'isPrintInfo', 1);
    parse(p, varargin{:});
    ip = p.Results;

    % 最大阶数
    N_MAX = max(n)+3;
    arg_row = arg(:).';

    if (abs(max(arg(:))) > 1.4e3) && (N_MAX>max(arg(:))/2)
        ip.Method = 'Recurrence';
    end

    switch ip.Method
        case 'MATLAB'
            J = besselj((0.5:1:N_MAX+0.5).', arg_row);
            j = sqrt(pi/2./arg_row) .* J;
            j = j(n+1,:);
        case 'Recurrence'
            N_MAX = max(n) + 2e2;
    %         N_MAX = max(n(:)) + round(2*abs(max(arg_row(:))))+9;
            j = ones(N_MAX+1, length(arg_row));
            j(N_MAX+1, :) = 1e-300;
            j(N_MAX, :) = 1e-300;
            for nn = N_MAX-1:-1:1
                j(nn,:) = j(nn+1,:).*(2*nn+1) ./ arg_row - j(nn+2,:);
            end
            j0 = sinc(arg_row/pi) ./ j(1,:);
            j = j .* j0;

            % 接近于零是，有特别解
            j(:, abs(arg_row)<1e-200) = 0;
            j(1, abs(arg_row)<1e-200) = 1;
            j = j(n(:)+1, :);
        case 'Normaliarged'
            jbar = SphBessel1Norm(n, arg_row);
    % 		j = jbar .* arg_row.^n .* sqrt(pi) ./ 2.^(n+1) ./ gamma(n+3/2);
            j = jbar .* sqrt(pi) .* exp(n.*log(arg_row) ...
                - (n+1).*log(2) - gammaln(n+3/2));
        otherwise
            error('错误的计算方法！\n');
    end

    j(n==0,arg_row==0) = 1;
    j(n~=0,arg_row==0) = 0;
    j = reshape(j, size(0*n.*arg));

end

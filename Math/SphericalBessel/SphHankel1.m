% ==========================================================================
% FUNCTION
%       计算球Hankel函数
% --------------------------------------------------------------------------
% 输入
% 	n	--- 阶数；必须在第1维
%	z	--- 宗量；必须在第2维及以上
% --------------------------------------------------------------------------
% 输出
%	h	---	球Hankel函数
% ==========================================================================
function h = SphericalHankel(n, zArgument, varargin)

    p = inputParser();
	addParameter(p, 'isPrintInfo', 1);
    parse(p, varargin{:});
    ip = p.Results;
    
    n = n(:);
    N = max(n);
    zArgument_row = zArgument(:).';

    if max(abs(imag(zArgument_row))) < 1e-200
        j = SphBessel1(n, zArgument);
        y = cal_sphNeumann(n, zArgument);
        h = j + 1i*y; % 有时不能用这种方法，因为j和y可能很相近，若相减会有很大误差
        return
    end

    h = zeros(N+1, length(zArgument_row));
    h(1,:) = exp(1i*zArgument_row)./(1i*zArgument_row);
    h(2,:) = exp(1i*zArgument_row).*(1-1i*zArgument_row) ./(1i*zArgument_row.^2);
    for nn = 1:N-1
        h(nn+1+1,:) = (2*nn+1)./zArgument_row .*h(nn+1,:) - ...
            h(nn-1+1,:);
    end
    h = h(n+1,:);
    h = reshape(h, size(0*n.*zArgument));

%     if ip.isPrintInfo
%         if check_infnan(h)
%             warning('出现inf或nan!\n');
%         end
%     end
end

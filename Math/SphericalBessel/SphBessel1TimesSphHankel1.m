%==========================================================================
% The spherical Bessel function of the first kind times 
%   the spherical Hankel function of the first kind.
% - j_n(zj) * h_n(zh)
%--------------------------------------------------------------------------
% INPUT
%   max_order
%           --- the maximum order
%   zj      --- the argument of the spherical Bessel function
%   zh      --- the argument of the spherical Hankel function
%--------------------------------------------------------------------------
% OUTPUT
%   res     --- The result, where the first dimension is the results 
%               corresponding to the degree 0, 1, 2, ..., max_order
%--------------------------------------------------------------------------
% dimension: order => zj .* zh
%==========================================================================

function res = SphBessel1TimesSphHankel1(max_order, zj, zh, varargin)
    p = inputParser;
    addParameter(p, 'method', 'robust');
    parse(p, varargin{:});
    ip = p.Results;
    
    zj_row = zj + 0 * zj .* zh;
    zh_row = zh + 0 * zj .* zh;
    zj_row = zj_row(:).';
    zh_row = zh_row(:).';

    switch ip.method
        case 'direct'
            n = (0:max_order).';
            j = SphBessel1(n, zj_row);
            h = SphericalHankel(n, zh_row);
            res = j .* h;
        case 'norm'
            n = (0:max_order).';
            j_bar = SphericalBesselJNorm(n, zj_row);
            h_bar = SphericalHankelNorm(n, zh_row);

%             res = 1 ./ (1i * (2*n+1) .* zh) ...
%                 .* exp(n .* log(zj ./ zh)) ...
%                 .* j_bar .* h_bar;
            res = 1 ./ (1i * (2*n+1) .* zh_row) ...
                .* exp(n .* log(zj_row ./ zh_row) + log(j_bar .* h_bar));
            
%              res(zj == 0 & n ~= 0 + 0*zj.*zh) = 0;
%             res(zj == 0 & n == 0 + 0*zj.*zh) = 1 ./ (1i .* zh) ...
%                 .* j_bar(n == 0) .* h_bar (n == 0);

        case 'robust'
            n = (0:max_order).';
            j_bar = SphBessel1Norm(n, zj_row);
            h_bar = SphHankel1Norm(n, zh_row);
            jh = SphBessel1(n, zj_row) .* SphHankel1(n, zh_row);
            res = 1 ./ (1i * (2*n+1) .* zh_row) ...
                .* exp(n .* log(zj_row ./ zh_row) + log(j_bar .* h_bar));
            
            zj_expand = zj_row + 0 * zj_row.*zh_row.*n;
            res(zj_expand == 0 & n ~= 0) = 0;
            res(zj_expand == 0 & n == 0) = 1 ./ (1i .* zh_row) ...
                .* j_bar(zj_row == 0 & n == 0) .* h_bar (n == 0);
            
            res(isnan(j_bar) | isinf(j_bar) | isnan(h_bar) | isinf(h_bar)) = ...
                jh(isnan(j_bar) | isinf(j_bar) | isnan(h_bar) | isinf(h_bar));
            
            res(:, abs(zj_row) > 1200 | abs(zh_row) > 1200) = ...
                jh(:, abs(zj_row) > 1200 | abs(zh_row) > 1200);
            
    end
    
    res = reshape(res, size(n .* zj .* zh));
end
